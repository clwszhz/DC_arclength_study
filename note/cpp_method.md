# C++ 求解方法简明笔记

这份笔记只记录当前 demo 的核心方法和源码对应关系，方便回看代码时定位。

## 1. 总体流程

- 程序入口在 `src/main.cpp`。
- 先解析网表，得到 `Circuit`。
- 如果没有 BJT，走线性 MNA：`A*x=b`。
- 如果含有 BJT，走非线性 MNA：`F(x)=0`。
- 线性和非线性最后都会输出变量名、解向量和残差。

```text
网表 -> Circuit -> 判断是否含 BJT
无 BJT: assembleMna() -> solveLup()
有 BJT: makeNonlinearMnaSystem() -> solvePseudoArclength()
```

## 2. 线性 MNA

- 代码位置：`src/Mna.cpp`。
- 输入：`Circuit`。
- 输出：`MnaSystem`，里面存 `A`、`b`、变量名和 stamp 日志。
- 未知量顺序：

```text
x = [非地节点电压, 电压源支路电流]
```

- 电阻 stamp：对角加 `+g`，非对角加 `-g`。
- 电流源 stamp：只改右端项。
- 电压源 stamp：增加一个支路电流未知量，并增加约束行。

## 3. LUP 求解

- 代码位置：`src/Lup.cpp`。
- 作用：求解线性方程 `A*x=b`。
- 使用行主元选取，避免主元太小导致数值不稳定。
- 线性电路直接调用一次 LUP。
- 非线性 Newton 和伪弧长 corrector 中也反复调用 LUP。

## 4. 非线性 MNA

- 代码位置：`src/Nonlinear.cpp`。
- 非线性方程不存成字符串表达式。
- 程序保存的是电路结构和变量顺序。
- 每次给定当前 `x`，现场计算：

```cpp
Vector evaluateResidual(circuit, system, x); // F(x)
Matrix evaluateJacobian(circuit, system, x); // J(x)
```

- `F(x)` 是 KCL/电压源约束/BJT 电流残差。
- `J(x)` 是 `F(x)` 对所有未知量的偏导矩阵。
- BJT 使用简化 Ebers-Moll 模型，计算 `Ic/Ib/Ie` 以及对应偏导。

## 5. 固定点同伦

- 代码位置：`src/Homotopy.cpp`。
- 目标：把难解的 `F(x)=0` 变成从简单问题逐步走向真实问题。
- 同伦方程：

```text
H(lambda,x) = (1-lambda)*gleak*(x-a) + lambda*Fh(x)
```

- `lambda=0` 时：

```text
H(0,x) = gleak*(x-a)
解为 x=a
```

- `lambda=1` 时：

```text
H(1,x) = Fh(x)
等价于原始电路方程 F(x)=0
```

- `a` 来自网表 `.START`，用于给同伦路径一个已知起点。
- `Fh(x)` 是行重排后的残差，用来对齐 MATLAB 参考代码的固定点同伦结构。

## 6. 伪弧长法

- 当前代码全程使用伪弧长，不再先尝试自然参数 continuation。
- 伪弧长把 `lambda` 也当成未知量：

```text
y = [lambda, x1, x2, ..., xn]
```

- 每一步分三件事：

```text
1. 求当前同伦曲线切线 tangent
2. predictor: y_pred = y + ds*tangent
3. corrector: Newton 修正到 H(lambda,x)=0
```

- corrector 解的是增广方程：

```text
H(lambda,x) = 0
tangent^T * (y - y_pred) = 0
```

- 第一行保证点在同伦曲线上。
- 第二行保证修正点不要跑离预测点太远。
- 如果 corrector 失败，弧长步长减半后重试。
- 如果 corrector 很快收敛，下一步步长适当放大。

## 7. 关键函数对应

- `evaluateHomotopy()`：计算 `H(lambda,x)`。
- `evaluateHomotopyJacobianX()`：计算 `dH/dx`。
- `evaluateHomotopyJacobianLambda()`：计算 `dH/dlambda`。
- `initialTangent()`：计算初始切线。
- `tangentAt()`：计算当前点切线，并保持方向连续。
- `makeAugmentedJacobian()`：生成 corrector 的增广 Newton 矩阵。
- `correctPseudoArclength()`：执行伪弧长 Newton corrector。
- `solvePseudoArclength()`：同伦主循环。
- `solveAtLambdaOne()`：到达 `lambda=1` 后，用原始 `F(x)=0` 做最终 Newton。

## 8. `solvePseudoArclength()` 主循环

主要逻辑：

```text
1. 检查 step_size、gleak 和 start 维度
2. 记录第 0 步：lambda=0, x=start
3. 生成 y=[lambda,x]
4. 计算初始切线
5. predictor 预测新点
6. corrector 修正新点
7. 成功则记录路径点
8. 判断当前路径段是否穿过 lambda=1
9. 如果穿过 lambda=1，插值得到交点，再做最终 Newton
10. 如果这个 DC 解不是重复解，就加入 result.solutions
11. 更新切线，继续追踪后面的路径
12. 超过最大步数或步长过小则停止
```

最核心的三句：

```cpp
y_predicted = addScaled(y, arc_step, tangent);
corrected = correctPseudoArclength(...);
tangent = tangentAt(...);
```

对应：

```text
预测 -> 修正 -> 更新切线
```

现在代码不会在第一次到达 `lambda=1` 时直接结束，而是继续沿同伦曲线走。每次路径从 `lambda<1` 到 `lambda>1`，或者反过来从 `lambda>1` 到 `lambda<1`，都说明同伦路径穿过了真实电路平面 `lambda=1`，程序就会记录一个候选 DC 解。

## 9. 为什么还要最终 Newton

- 伪弧长路径穿过 `lambda=1` 时，插值得到的点只是接近真实电路解。
- 程序最后固定 `lambda=1`，直接解原始方程：

```text
F(x)=0
```

- 这一步用标准 `evaluateResidual()` 和 `evaluateJacobian()`，不是行重排后的 `Fh/Jh`。
- 最终输出的 `max|F|` 才是判断真实电路工作点是否成立的关键。
- 如果最终 Newton 不收敛，或者新解和已有解太接近，程序不会把它加入解列表。

## 10. 当前验证

默认参数：

```text
step_size = 0.1
gleak = 1e-3
tolerance = 1e-9
max_steps = 3000
```

验证结果：

```text
schmitt1.cir: 找到 3 个 lambda=1 DC 解，max|F| 约 1e-11
schmitt2.cir: 找到 3 个 lambda=1 DC 解，max|F| 约 1e-11
chua.cir:     找到 9 个 lambda=1 DC 解，max|F| 约 1e-11
ctest: 100% passed
```

这里的“多个解”来自 `HomotopyResult::solutions`。每个元素保存最终 Newton 后的 `x`、残差范数和 Newton 迭代次数。

## 11. 当前不足

- 还没有稀疏矩阵。
- 还没有 AC、TRAN、DC sweep。
- BJT 模型只覆盖当前示例需要的参数。
- 伪弧长是教学版，没有完整 MATLAB `pchomotopy` 的高阶可变步长。
- 当前“全解”指沿当前固定点同伦路径穿过 `lambda=1` 的所有已记录交点，不等于数学上证明所有孤立分支都已经找到。
