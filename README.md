# DC 工作点求解 demo

这是一个用于学习电路 DC 工作点求解的 C++17 demo。核心目标是把网表自动转成 MNA 方程，再用 LUP、Newton 和伪弧长同伦求解 DC 工作点。

## 当前功能

- 线性 `.OP`：支持 `R/I/V`，组装 `A*x=b`，用手写稠密 LUP 求解。
- 非线性 `.OP`：支持 `R/I/V/Q`，其中 `Q` 是简化 NPN Ebers-Moll BJT。
- 支持最小网表语法：`R`、`I`、`V`、`Q`、`.MODEL`、`.START`、`.OP`、`.END`。
- 非线性电路不手写某个 case 的公式，而是根据网表自动计算 `F(x)` 和 `J(x)`。
- 伪弧长路径每次穿过 `lambda=1` 时都会做最终 Newton，并把不同 DC 解记录下来。

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
     -> 记录所有穿过 lambda=1 的 DC 解
```

固定点同伦形式：

```text
H(lambda,x) = (1-lambda)*gleak*(x-a) + lambda*Fh(x)
```

`a` 来自 `.START`。`Fh(x)` 是为了对齐 MATLAB 参考程序而做过行重排的 MNA 残差。伪弧长法把 `y=[lambda,x]` 当成整体变量，先沿切线预测，再用 Newton corrector 解：

```text
H(lambda,x) = 0
tangent^T * (y - y_pred) = 0
```

这样即使 `lambda` 发生折返，也可以继续沿同伦曲线追踪。

## 编译运行

示例使用 MSYS2 UCRT64：

```bash
cd /d/VsCodeP/study/new_reading/demo
cmake -S . -B build-msys -G Ninja
cmake --build build-msys
```

运行示例：

```bash
./build-msys/dcsolve.exe examples/divider.cir
./build-msys/dcsolve.exe examples/schmitt1.cir
./build-msys/dcsolve.exe examples/schmitt2.cir
./build-msys/dcsolve.exe examples/chua.cir
```

运行测试：

```bash
ctest --test-dir build-msys --output-on-failure
```

## 当前验证结果

默认参数为 `step_size=0.1`、`gleak=1e-3`、`max_steps=3000`。当前运行结果：

```text
schmitt1.cir: 找到 3 个 lambda=1 DC 解，max|F| 约 1e-11
schmitt2.cir: 找到 3 个 lambda=1 DC 解，max|F| 约 1e-11
chua.cir:     找到 9 个 lambda=1 DC 解，max|F| 约 1e-11
dcsolve_tests: passed
```

## 目录说明

```text
include/dcsolve/   头文件和数据结构
src/               Parser、MNA、Nonlinear、Homotopy、LUP 和 main
examples/          可运行网表示例
tests/             LUP 单元测试
note/              学习笔记
```

## 限制

- 目前只做 DC 工作点，不做 AC、TRAN、DC sweep。
- 矩阵使用稠密存储，还没有稀疏矩阵。
- BJT 模型只覆盖示例需要的简化参数。
- 当前“多解”是指沿这条固定点同伦路径所有穿过 `lambda=1` 的交点；这不是对所有可能孤立分支的数学证明。
