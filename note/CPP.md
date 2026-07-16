# C++ 语法与工程笔记

## 1. 头文件 `.h` 与源文件 `.cpp` 的分工

在 C++ 工程中，`.h` 头文件主要写声明，说明一个模块提供了哪些类型和函数；`.cpp` 源文件主要写实现，说明这些函数具体怎么执行。其他文件想使用某个模块时，通常 `#include` 这个模块的 `.h`，而不是包含它的 `.cpp`。

`main.cpp` 一般不需要 `main.h`，因为 `main.cpp` 通常是程序入口，负责调用其他模块并组织主流程；只有当某些函数或类型需要被其他文件复用时，才有必要写到头文件里。

```text
例子：
Circuit.h    声明 Circuit、Resistor 等数据结构
Circuit.cpp  实现 Circuit 相关函数
main.cpp     #include "dcsolve/Circuit.h" 后使用这些声明
```

## 2. 命名空间 `namespace`

`namespace` 用来把名字分组，避免不同模块中的函数、变量或类型重名。即使已经 `#include` 了头文件，如果函数或类型声明在某个命名空间中，使用时仍然要加上命名空间名。

匿名命名空间 `namespace { ... }` 没有名字，里面的函数或变量只在当前 `.cpp` 文件中可见，常用于放当前文件内部使用的小工具函数。

两个命名空间可以嵌套使用，例如 `namespace dcsolve { namespace { ... } }`：外层表示属于 `dcsolve` 模块，内层匿名命名空间表示这些辅助函数只在当前 `.cpp` 文件内部可见。

```cpp
namespace dcsolve {
    void solve();
}

namespace {
    void helper();
}

int main() {
    dcsolve::solve(); // 使用 dcsolve 命名空间中的函数
    helper();         // 只能在当前 cpp 文件中使用
}
```

## 3. `if` 与 `try-catch`

`if` 用于条件判断：当条件成立时执行对应代码，否则跳过。在 `main.cpp` 中常用于检查输入是否合法，例如参数个数不对就打印用法并退出。

`try-catch` 用于异常处理：把可能出错的核心流程放进 `try`，如果其中抛出异常，就由 `catch` 接住并统一处理错误。

```cpp
if (argc != 2) {
    printUsage(argv[0]);
    return 2;
}

try {
    auto circuit = dcsolve::parseNetlistFile(argv[1]);
} catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return 1;
}
```

## 4. 命令行参数与 `main(argc, argv)`

在终端中输入命令时，格式通常是“程序名 + 参数”。操作系统启动可执行程序后，程序从 `main` 函数开始执行，命令行中的内容会传给 `main(int argc, char** argv)`。

`argc` 表示参数个数，`argv` 保存每个参数字符串。`argv[0]` 通常是程序名，后面的 `argv[1]`、`argv[2]` 才是用户传入的参数。

```text
命令：
./dcsolve.exe examples/divider.cir

对应：
argc = 2
argv[0] = "./dcsolve.exe"
argv[1] = "examples/divider.cir"
```

## 5. `std::ifstream` 与 `std::getline`

`std::ifstream` 是文件输入流，用来从文件中读取数据。`std::ifstream input(path);` 表示打开 `path` 指向的文件，并用 `input` 作为读取文件的通道。

`std::getline(input, line)` 从输入流中读取一整行文本，存入 `line`。读取成功返回 `true`，读到文件末尾或失败返回 `false`，所以常用于 `while` 循环逐行处理文件。

```cpp
std::ifstream input(path);
std::string line;

while (std::getline(input, line)) {
    // 每次循环处理文件中的一行
}
```

## 6. `std::string::find` 与 `substr`

`find` 用来在字符串中查找某个字符或子串，找到就返回第一次出现的位置；没找到就返回 `std::string::npos`。常用它判断一行文本中是否含有某个标记。

`substr(pos, len)` 用来截取字符串的一部分。在网表解析中，可以用 `find(';')` 找到注释开始位置，再用 `substr(0, pos)` 保留注释前面的有效内容。

```cpp
auto pos = line.find(';');
if (pos != std::string::npos) {
    line = line.substr(0, pos);
}
```

## 7. `throw` 抛出异常

`throw` 用来抛出异常，表示当前代码遇到了错误，不能继续正常执行。它会立刻中断当前函数，并把错误向外层传递，直到被某个 `catch` 接住。

在本工程中，解析网表出错时会 `throw std::runtime_error(...)`，最后由 `main.cpp` 里的 `catch` 统一打印错误并退出程序。

```cpp
if (words.size() != expected) {
    throw std::runtime_error("wrong field count");
}
```

## 8. `std::vector<T>`

`std::vector<T>` 是动态数组，`T` 表示数组中存放的元素类型。它可以保存多个同类型对象，并且大小可以随 `push_back` 等操作增长。

在本工程中，`std::vector<Resistor> resistors;` 表示保存多个电阻对象，解析网表时每读到一个电阻，就可以把它加入这个列表。

```cpp
std::vector<int> nodes;
std::vector<Resistor> resistors;
resistors.push_back(r);
```

## 9. 成员函数末尾的 `const`

成员函数参数列表后面的 `const` 表示该函数不会修改当前对象。它常用于只读取对象内容、但不改变对象状态的函数。

例如 `std::vector<int> Circuit::nodeLabels() const` 表示 `nodeLabels` 是 `Circuit` 的只读成员函数，可以读取电路中的节点信息，但不能修改 `Circuit` 内部的元件列表。

```cpp
std::vector<int> Circuit::nodeLabels() const {
    // 只读取 Circuit，不修改 Circuit
}
```

## 10. `std::set`、范围 `for` 与 lambda

`std::set<T>` 是集合，会自动去重并排序。范围 `for` 可以直接遍历容器中的每个元素。lambda 是临时定义的小函数，常用于当前函数内部的简单操作。

在 `Circuit::nodeLabels()` 中，`set<int>` 用来收集所有非地节点并去重；`add_node` 是 lambda，用来判断节点是否为 0；多个范围 `for` 分别遍历电阻、电流源和电压源。

```cpp
std::set<int> labels;
auto add_node = [&labels](int node) {
    if (node != 0) labels.insert(node);
};

for (const auto& r : resistors) {
    add_node(r.a);
    add_node(r.b);
}
```

## 11. 成员访问、`reserve` 与 `assign`

连续的点号表示逐层访问成员，例如 `system.variable_names.reserve(n)` 是先访问 `system` 的 `variable_names`，再调用这个 `vector` 的 `reserve` 函数。`reserve(n)` 用来提前预留空间，减少后续插入时的扩容。

`assign(count, value)` 会清空容器原内容，并重新放入 `count` 个 `value`。在 MNA 初始化中，常用它生成全 0 的矩阵和向量。

```cpp
system.variable_names.reserve(n);
system.b.assign(n, 0.0);
system.a.assign(n, Vector(n, 0.0));
```
