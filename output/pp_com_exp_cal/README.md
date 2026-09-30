# pp 弹性散射九能量实验与 Fortran 对比

能量为 **3.037、9.69、14.16、50.06、98、142、239.9、270.8、286.7 MeV**（入射质子实验室动能）；共有190个EXFOR实验点。质心系角度和微分截面采用原记录定义，不镜像复制数据、不施加相同粒子计数因子、不拟合理论归一化。

## 文件

- `pp_experimental_angular_3x3.png/.pdf`：实验数据3×3图，PNG为高分辨率，PDF为矢量。
- `pp_exp_cal_angular_3x3.png/.pdf`：实际 Fortran 输出与实验对比图（已完成）。
- `experimental_points.csv`、`pp_exp_*MeV.csv`：全部/分能量作图数据，保留原单位、原值、误差、EXFOR编号和来源行号。
- `queries/`：九份原始请求、QueryOutcome、原始展开表和标准化prepared_panel。没有直接解析或修改原始`.x4`记录。
- `plot_pp_3x3.py`：由冻结查询结果重画实验图。
- `plot_pp_comparison.py`：由冻结查询结果和 Fortran 输出重画比较图。
- `query_experiments.py`：需要本地EXFOR数据库；重新运行九份QuerySpec。
- `exfor_tool/`：绘图/查询工具源码快照；正式源码在工作区 `tools/EXFOR`。
- `panel_queries.json`：能量与原始查询输出清单。
- `DATA_SOURCES.md`：每组论文、角度覆盖、误差与筛选规则；使用前应阅读。
- `plot_manifest.json`、`comparison_manifest.json`：实验/比较图的完整数据和处理记录；比较图记录实际运行 `pp_summary.txt` 参数。
- `tool_provenance.json`、`data_validation.json`：工具源码哈希、数据库版本和数据验证记录。

## 重画

在此目录中执行（Python仅负责画图，物理解算由Fortran完成）：

```bash
python3 -m pip install -r requirements.txt
python3 plot_pp_3x3.py
python3 plot_pp_comparison.py --batch fortran_runs/实际批次时间
```

比较脚本读取每个批次下 `<能量>MeV/pp_angular_distribution.csv`（例如`3p037MeV`）与 `pp_summary.txt`，CSV字段为`theta_cm_deg`和`dsigma_domega_mb_sr`。不指定`--batch`时读取`latest_fortran_batch.txt`。具体计算设置、验证与局限另见Fortran计算说明；本脚本不硬编码未验证的网格或库仑半径。

图横轴0°–180°，逐图独立纵轴。库仑前后向尖峰会发散；比较图的纵轴范围依据实验值及理论在实测角域的数值设定，超出显示上限的尖峰会裁切。原始理论CSV仍保留全部输出角度，不能把图的有限显示范围当成总截面截止。

实验误差不全是1σ独立误差，多组还带共同尺度误差；本图用于形状和绝对量级比较，不据此作未考虑相关性的χ²结论。

计算设置、数值验证、独立参考及误差边界见[FORTRAN_COMPARISON.md](FORTRAN_COMPARISON.md)。

纵轴已按用户要求改为对数坐标（截面数值及单位mb/sr不变）；各子图范围独立，PNG/PDF与重画代码同步更新。

## 可选：重新查询EXFOR

本目录已保留当前图使用的查询结果与实验数表。运行 Fortran 和离线重画均不需要 EXFOR 数据库，也不需要原来的工作区路径。只有希望重新访问数据库时才使用 `query_experiments.py`；重新查询会更新 `queries/` 中的数据，应先保存旧版本。

数据库和SQLite索引单独安装，以明确路径传入：

```sh
python3 query_experiments.py \
  --index /path/to/exfor.sqlite \
  --source-root /path/to/exforall
```

默认调用本目录的 `exfor_tool/` 源码快照。如需使用另外安装的工具实现，可再加 `--tool-src /path/to/EXFOR/src`（该目录必须含 `exfor_tool/query.py`）；兼容原工作区布局的简写为 `--workspace /path/to/workspace`。来源JSON中保存的原机器绝对路径仅用于记录历史查询来源，不参与离线重画。
