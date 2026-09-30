# 氘核：实验与计算比较

首先打开 `氘核_实验与计算比较.pdf`。第1页是比较表，第2页为密度图与半径定义，第3页给出来源和复现信息。

- `comparison_table.csv`：6项计算—实验比较。点核子半径移至正文说明；四极矩、磁矩采用单核子冲量近似；电荷半径是有限核子尺寸加Darwin–Foldy项的估计。
- `deuteron_density.png/pdf`：相对距离径向概率与未极化、球平均点核子总数密度。密度归一为2，不是实验电荷密度。
- `calculation_values.json`：从Fortran输出读取的数值，含所用附加电磁常数。
- `input.nml`、`run_deuteron.sh`：本次192点动量网格、60 fm坐标上限、0.01 fm步长的输入与计算脚本。
- `latest_run.txt`、`fortran_runs/`：最新实际运行位置以及逐次运行输出；内含 `deuteron_observables.dat` 和 `deuteron_coordinate_wavefunction.csv`。
- `experimental_sources/`：实验/评价数值、性质、原始来源与历史下载校验记录；第三方PDF获取地址见 [SOURCE_LINKS.md](SOURCE_LINKS.md#papers)。
- `validation/`：数值收敛和独立半径核验。更完整开发审计位于项目 `runs/deuteron_observables_2026-09-29/`。
- `build_report.py`、`氘核_实验与计算比较.tex`：图和PDF生成源；`build/`保存排版过程文件。
- `provenance.json`、`SHA256SUMS`：交付时的计算来源和文件校验记录。

在本目录运行：

```bash
./run_deuteron.sh             # 调用 ../../app/fewbody.x，必要时先make
python3 build_report.py      # 读取最新Fortran结果，生成图和PDF
```

计算需要GNU Fortran和Make，不依赖Python；图表生成需要Python、Matplotlib和XeLaTeX/ctex。程序源码位于 `../../app/deuteron_observables.f90`，算法说明见[《EM500核力与两核子态求解》](../../docs/EM500核力与两核子态求解.pdf)第3章（第3.4—3.5节：坐标波函数、半径与其他观测量）。统一代码根目录也可运行 `./run.sh input_deuteron_observables.nml`，自动按结束时间归档新输出。

现有EM500及求解器未重拟合。现代结合能评价与历史定参目标不同；完整电荷/磁形状因子及双核子电磁算符尚未实现，不能以本表的近似电磁量直接宣称精密核力验证已完成。改动输入并重算后，旧校验清单不再对应新输出，应另存新的研究记录。
