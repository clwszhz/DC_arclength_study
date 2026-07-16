# MATLAB 阅读笔记

## 1. MATLAB 的下标从 1 开始

MATLAB 的数组和矩阵下标从 `1` 开始，而不是从 `0` 开始。因此 MATLAB 中的

```matlab
fc(1,1)
```

对应 C/C++ 或 Python 这类语言中的

```text
fc[0][0]
```

如果 `fc` 是列向量，`fc(1,1)` 就是第一行第一列的元素，也可以简写成 `fc(1)`。

```matlab
fc = 0 * y0;
fc(1,1) = 1;
```

这表示把列向量 `fc` 的第一个元素设成 `1`。

## 2. `abs`、`|`、`&` 与逻辑判断

`abs(x)` 表示取绝对值，常用于判断一个数是否足够接近目标值。`|` 表示逻辑“或”，左右条件至少一个为真则结果为真；`&` 表示逻辑“与”，左右条件都为真时结果才为真。

在同伦路径追踪中，`close` 判断当前点是否接近 `lambda = 1`，`cross` 判断一步是否跨过 `lambda = 1`。其中 `tol(2)` 表示容差向量 `tol` 的第 2 个元素。

```matlab
close = (abs(lambdatemp - 1) < tol(2)) | (abs(lambda - 1) < tol(2));
cross = ((lambda < 1) & (lambdatemp > 1)) | ...
        ((lambdatemp < 1) & (lambda > 1));
```
