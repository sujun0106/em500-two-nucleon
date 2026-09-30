# np角分布：九能量实验与Fortran EM500对比

当前能量为3.1、10.04、25.8、50、96、137、150、200、300.2 MeV。按用户要求，原284 MeV被300.2 MeV替换，上限相应放宽到300.2 MeV。全图共242个实验点。

## 图与数据

- `np_exp_cal_angular_3x3.png` / `.pdf`：实验—Fortran计算对比图。
- `np_experimental_angular_3x3.png` / `.pdf`：纯实验图。
- `experimental_points.csv`、`np_exp_*MeV.csv`：实际绘图数据，含逐点误差、来源及原始/有效单位。
- `fortran_runs/`：实际Fortran运行结果；`latest_fortran_batch.txt`指向当前图所用批次。新批次前八能量沿用上一批结果，300.2 MeV新算721角度；`assembly_provenance.json`记录复用来源。
- `fortran_inputs/`、`run_fortran_nine.sh`：九能量输入与重新计算入口。
- `comparison_at_experimental_angles.csv`：计算曲线线性插值到实验角度的差值表，未作实验拟合。
- `FORTRAN_COMPARISON.md`、`fortran_summary.csv`、`fortran_provenance.json`：计算设置、数值检查和源码哈希。

## 实验来源

- 3.1 MeV：6点；20319002。
- 10.04 MeV：6点；13782002。
- 25.8 MeV：16点；10355002, 10355003。
- 50 MeV：20点；10355002, 10355003。
- 96 MeV：32点；22223002。
- 137 MeV：14点；21800002, 21800003。
- 150 MeV：16点；11123003。
- 200 MeV：21点；41209002。
- 300.2 MeV：111点；Franz 2000。

前八能量来自EXFOR，数据commit为`211d7b5195fa426b40d29393b642975f75a0542a`。查询及原始数表保存在queries/；数字格式解析修正和历史筛选记录已归档在项目runs/中。

300.2 MeV来自Franz等，Physica Scripta T87 (2000) 14–21，DOI: https://doi.org/10.1238/Physica.Topical.087a00014；NN-OnLine组2250–2253。111点覆盖79.8–179.2度，全部保留，包含子组重叠角度。并非EXFOR记录：当前精确能量查询仍返回no_match，响应文件保留在queries/。

NN-OnLine相对数表乘论文Table I给出的10.63 mb/sr，得到已冻结的绝对表。对应逐点误差同比例换算；论文给出的3.1%共同归一化误差单独保留，不能当作111个独立误差。该绝对化有论文依据，不是为贴合理论曲线而调整单位或归一化。来源说明、相对表和派生绝对表保存在franz_300p2/；历史原文校验信息随包保留，论文与网页原件从[SOURCE_LINKS.md](SOURCE_LINKS.md#franz)所列原站获取。CSV中Franz行的data_original_value指已冻结的派生绝对表值，原相对值在franz_300p2/tables/FR00_g*.csv。

Franz实验使用20 MeV能量分箱，本图理论为300.2 MeV名义能量计算，尚未按真实束流能谱卷积。所有角度为质心系出射中子角度，截面mb/sr。每幅图纵轴刻度独立，误差棒不包含共同归一化误差。

## 单位修正与历史版本

当前unit_overrides.json为空；284 MeV数据、经验单位修正规则及旧图已归档到项目runs/np_comparison_replace300p2_2026-09-29/before/。当前图不再采用该修正。其他八个能量的实验数据和计算数值不变。

## 复现

在本目录运行（Python需Matplotlib）：

```sh
python3 plot_np_3x3.py
python3 plot_np_comparison.py
```

仅重画无需EXFOR数据库或重复Fortran计算。重新运行九能量计算可执行`./run_fortran_nine.sh`，之后重画对比图。原Fortran程序没有改动。

`query_experiments.py`可重新查询本地EXFOR；对300.2 MeV保留no_match，同时继续使用明确标识的Franz来源，不伪造EXFOR记录。`panel_queries.json`的prepared_panels字段指定外部实验来源；`franz_300p2/prepared_panel.json`保留绘图记录和来源说明。

## 可选：重新查询EXFOR

本目录已保留当前图使用的查询结果与实验数表。运行 Fortran 和离线重画均不需要 EXFOR 数据库，也不需要原来的工作区路径。只有希望重新访问数据库时才使用 `query_experiments.py`；重新查询会更新 `queries/` 中的数据，应先保存旧版本。

数据库和SQLite索引单独安装，以明确路径传入：

```sh
python3 query_experiments.py \
  --index /path/to/exfor.sqlite \
  --source-root /path/to/exforall
```

默认调用本目录的 `exfor_tool/` 源码快照。如需使用另外安装的工具实现，可再加 `--tool-src /path/to/EXFOR/src`（该目录必须含 `exfor_tool/query.py`）；兼容原工作区布局的简写为 `--workspace /path/to/workspace`。来源JSON中保存的原机器绝对路径仅用于记录历史查询来源，不参与离线重画。
