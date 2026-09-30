# 随包提供的 DSYEV / DGESV 数值求解器

本目录是官方 Reference-LAPACK **v3.12.1** 中 `DSYEV` 和 `DGESV` 的源码依赖闭包，包含 60 个未修改的 Fortran 源文件。随包提供对应参考 BLAS 子程序，因此默认构建不需要服务器预装 BLAS/LAPACK/MKL，也不需要 Python、C 或 C++ 编译器。编译器正常自带的 Fortran 与系统运行库仍需可用。

- 官方发布页：https://github.com/Reference-LAPACK/lapack/releases/tag/v3.12.1
- 仓库：https://github.com/Reference-LAPACK/lapack
- 固定 commit：`6ec7f2bc4ecf4c4a93496aa2fa519575bc0e39ca`
- 原始许可：BSD-3-Clause，见 `LICENSE`。
- 上游文件路径、固定 commit 下载 URL 和逐文件 SHA-256：见 `UPSTREAM.json`。
- 检查上游文件完整性：在本目录运行 `sha256sum -c SHA256SUMS`。
- `README.md`、`sources.mk`、`UPSTREAM.json`、`SHA256SUMS` 为本项目新增的封装/核验文件；开发测试代码移到根目录的 `other/tests/`；上游计算例程没有修改。

## 编译与链接

氘核本征问题和散射 K 矩阵本征相移使用双精度实对称本征求解接口 `DSYEV`；np 散射主值离散方程使用带部分选主元 LU 分解的 `DGESV`。全部依赖文件列在 `sources.mk` 的 `LAPACK_SOURCES` 变量中，路径相对于 `code/`。

建议 GNU Fortran 选项为 `-O2 -ffp-contract=off -frecursive`。不要使用 `-ffast-math`，不要使用 `-fdefault-real-8` 或 `-fdefault-integer-8`。本包采用正常的双精度实数和默认 32 位 LAPACK 整数接口；若将来链接其他数值库，必须使用一致的整数接口。

模块依赖顺序如下；并行 `make` 必须显式表达这些依赖：

1. `SRC/la_constants.f90`。
2. `SRC/la_xisnan.F90` 依赖 `la_constants`；`SRC/dlartg.f90` 依赖 `la_constants`。
3. `SRC/dlassq.f90` 依赖 `la_constants` 和 `la_xisnan`。
4. 其他源文件没有 Fortran `USE` 模块依赖。

编译时用 `-J` 指定模块输出目录，并用 `-I` 将该目录加入搜索路径。扩展名为 `.F90` 的文件须保留大写 F，以使 GNU Fortran 正常预处理。无需启用 OpenMP。

这些对象可直接加入最终链接命令，或用 `ar` 打包为静态库。请勿同时链接这套例程和另一份同名 BLAS/LAPACK 实现。

参考 BLAS 以可移植、可核验为目标，没有针对具体 CPU 调优。这里的氘核矩阵规模较小，默认优先采用随包源码；今后可选链接经过核验的服务器数值库。

## 独立数值检查

开发目录`other/tests/test_dsyev.f90`（路径相对code根目录） 测试 80 阶对称三对角矩阵：主对角为 2、副对角为 -1，其解析本征值为

`lambda_k = 2 - 2*cos(k*pi/81)`, `k = 1,...,80`。

测试覆盖上/下三角输入、只求本征值/同时求本征向量、总体缩放 `1`、`1e-180`、`1e180`，共 12 个组合；分别检查本征值、归一化本征向量残差和正交性。该测试仅依赖随包对象，无外部数值库。项目首次集成时使用 GNU Fortran 11.4.0，12 项均通过，各相应最大绝对误差不超过 `2.7e-15`（极端缩放的误差先除以总体缩放因子）。

np 散射扩展增加了同版本的 `DGESV`、`DGETRF`、`DGETRF2`、`DGETRS`、`DLASWP`、`IDAMAX`、`DTRSM`，其余依赖与已有 DSYEV 共用。`other/tests/check_np_solver.f90` 用有限积分上限下的可解析 Yamaguchi 势、常数势和旋转耦合双通道检验线性求解与主值消奇，避免仅以两份离散程序相符作为正确性依据。
