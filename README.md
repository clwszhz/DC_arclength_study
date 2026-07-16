# DC 工作点求解 demo

这是一个用于学习电路 DC 工作点求解的 C++17 demo。

## 当前功能

- 支持线性 `.OP`：`R/I/V` 元件，直接组装 `A*x=b`，用手写稠密 LUP 求解。
- 支持非线性 `.OP`：`R/I/V/Q`，其中 `Q` 为简化 NPN Ebers-Moll BJT。
- 支持最小网表语法：`R`、`I`、`V`、`Q`、`.MODEL`、`.START`、`.OP`、`.END`。
- 非线性电路不手写某个电路方程，而是根据网表自动生成 `F(x)` 和 `J(x)`。

## 核心方法

线性电路：

```text
网表 -> Circuit -> MNA stamp -> A*x=b -> LUP -> 节点电压/电压源电流
```

非线性电路：

```text
网表 -> Circuit -> 非线性 MNA -> F(x), J(x)
     -> 固定点同伦 H(lambda,x)
     -> 伪弧长 predictor-corrector
     -> lambda=1 后最终 Newton
```

不动点同伦形式为：

```text
H(lambda,x) = (1-lambda)*gleak*(x-a) + lambda*Fh(x)
```

其中 `a` 是 `.START` 给出的初始点，`Fh(x)` 是为对齐 MATLAB 参考程序而做过行重排的 MNA 残差。伪弧长法把

```text
y = [lambda, x]
```

作为整体变量，沿同伦曲线切线预测，再用 Newton corrector 解增广方程：

```text
H(lambda,x) = 0
tangent^T * (y - y_pred) = 0
```

这样不要求 `lambda` 每一步单调增加，路径发生折返时也能继续追踪。

## 编译运行

示例使用 MSYS2 UCRT64：

```bash
cd /d/VsCodeP/study/new_reading/demo
cmake -S . -B build-msys -G Ninja
cmake --build build-msys
```

运行线性例子：

```bash
./build-msys/dcsolve.exe examples/divider.cir
```

运行非线性 Schmitt 例子：

```bash
./build-msys/dcsolve.exe examples/schmitt1.cir
./build-msys/dcsolve.exe examples/schmitt2.cir
```

运行测试：

```bash
ctest --test-dir build-msys --output-on-failure
```

## 当前验证结果

默认参数下：

```text
schmitt1.cir: max|F| = 1.663995330e-12
schmitt2.cir: max|F| = 3.796460542e-12
dcsolve_tests: passed
```

## 目录说明

```text
include/dcsolve/   头文件和数据结构
src/               Parser、MNA、Nonlinear、Homotopy、LUP 和 main
examples/          可运行网表示例
tests/             LUP 单元测试
note/          学习笔记，不属于核心工程
```

## 限制

- 目前只做 DC 工作点，不做 AC、TRAN、DC sweep。
- 矩阵使用稠密存储，还没有稀疏矩阵。
- BJT 模型只覆盖老师示例需要的简化参数。
- 伪弧长实现是教学版 predictor-corrector，没有完整复现 MATLAB 高阶可变步长和多解记录。
