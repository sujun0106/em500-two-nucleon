# EM500核力与两核子态求解

Fortran 程序包 **v0.5.2**。使用固定的 EM500（Entem–Machleidt 2003 Idaho N³LO，500 MeV）两核子相互作用，计算氘核束缚态、np 弹性散射和 pp 弹性散射。程序还提供氘核波函数、点核子密度，以及指定算符近似下的半径和电磁观测量。

[完整说明书与三个操作范例](docs/EM500核力与两核子态求解.pdf) · [可编辑正文](docs/EM500核力与两核子态求解.md) · [版本记录](CHANGELOG.md) · [来源与许可](THIRD_PARTY_NOTICES.md)

## 获取、编译与运行

编译和物理计算使用 **GNU Fortran、GNU Make、Bash 及 GNU 基础命令**。建议在 Linux 或 Windows 的 WSL 中运行。首次准备 COUL90 库仑函数时，还需联网并安装 `curl` 或 `wget`。EM500 和所需的 LAPACK/BLAS 源码随仓库提供。

```sh
git clone https://github.com/sujun0106/em500-two-nucleon.git
cd em500-two-nucleon
make
make verify-sources
./run.sh input.nml
```

首次 `make` 会调用 `fetch_external.sh`，从原站获取固定版本的 COUL90，并校验 SHA256 后放入 `app/coulomb/coul90.f`。后续编译使用该本地文件。离线机器可从已有研究环境或原站取得文件，再执行：

```sh
./fetch_external.sh --source /path/to/COUL90.FOR
make
```

COUL90 的来源、作者和固定校验值见 [依赖说明](app/coulomb/README.md)。该文件由使用者从原来源获取，许可按原来源处理。

根目录提供三个体系的输入范例：

```sh
./run.sh input.nml                     # 氘核结合能与动量波函数
./run.sh input_deuteron_observables.nml # 氘核坐标波函数、半径与密度
./run.sh input_np.nml                  # np 弹性散射，96 MeV
./run.sh input_pp.nml                  # pp 弹性散射，50.06 MeV
```

计算先在根目录写结果，成功后连同实际输入和日志归档到 `output/YYYYMMDD_HHMMSS_UTC/`。同一根目录依次运行各任务；并行任务使用独立副本。第7章给出编辑输入、查看结果和归档的完整步骤。

统一可执行入口为 `app/fewbody.x`。输入采用 Fortran namelist：`calculation` 选择任务，`energy_lab_MeV` 指定入射粒子轰击静止靶核子的实验室动能。散射输出角度采用质心系，微分截面单位为 mb/sr。

## 主要输出与运行管理

| 任务 | 文件 | 内容 |
|---|---|---|
| 氘核 | `deuteron_summary.txt` | 本征能与结合能 `binding_MeV` |
| 氘核 | `deuteron_wavefunction.csv` | 动量空间 S/D 波函数 |
| 氘核可选 | `deuteron_observables.dat`、`deuteron_coordinate_wavefunction.csv` | 半径、单体电磁量、渐近参数、坐标波函数和点核子密度 |
| np | `np_angular_distribution.csv`、`np_summary.txt` | 角分布、输入参数与数值诊断 |
| pp | `pp_angular_distribution.csv`、`pp_summary.txt` | 角分布、输入参数与数值诊断 |

散射还保存分波、矩阵、振幅和总角动量截断检查。每次归档记录输入、版本、编译器与可执行文件信息、日志及校验清单。

```sh
./submit.sh input_np.nml          # 在当前机器后台运行
./run.sh --no-archive input_pp.nml # 完成后在根目录查看结果
./archive_run.sh                 # 手动归档
./remove.sh                      # 清理根目录的运行文件
./remake.sh                      # 清理 app 中的编译产物
```

`remove.sh` 会删除根目录尚未归档的已知结果；需要保留的任务请先归档。清理和重新编译应在当前任务结束后进行。`submit.sh` 使用 nohup；Slurm/PBS 服务器可在队列脚本中调用 `./run.sh 输入文件`。关闭计算机会中断任务，当前两体程序重启后从头运行。

默认编译选项与测试配套使用。数值接口依赖浮点精度和整数类型，修改 `-ffast-math`、`-fdefault-real-8`、`-fdefault-integer-8` 等选项需要重新验证。

## 文件结构

```text
app/                     Fortran 源码及本地编译产物
  em500/                 固定 EM500 实现与许可
  lapack/                线性代数源码与许可
  coulomb/               COUL90 来源、校验及本地下载位置
  tests/                 六个 Fortran 验证程序
    reference/           冻结参考数值和来源说明
  provenance/            核力配置、上游版本及校验记录
docs/                    主说明书及正文、TeX、排版脚本和素材
output/
  d_com_exp_cal/          氘核实验与计算比较
  np_com_exp_cal/         np 九能量实验与计算比较
  pp_com_exp_cal/         pp 九能量实验与计算比较
fetch_external.sh        获取并校验 COUL90
Makefile、input*.nml、*.sh 编译、输入、运行与清理入口
```

三个实验比较目录保存数值数据、来源记录、Fortran 输入与结果、自制图和绘图脚本。第三方论文全文及原网页由相应的 `SOURCE_LINKS.md` 提供获取地址。编译产物、下载的 COUL90、缓存及新计算的时间目录由 `.gitignore` 管理，留在本地。

## 验证与复现

在当前任务结束后，顺序执行：

```sh
make verify-sources
make check
```

第一条核验固定核力、数学依赖及参考表；第二条运行六项 Fortran 检查，覆盖氘核回归与解析观测量、np 解析模型与自旋振幅、pp 库仑匹配与交换处理。测试所需参考表位于 `app/tests/reference/`。编译与测试共用 `app/` 模块文件，按顺序执行。

根目录 `SHA256SUMS` 对应公开发布的文件清单。下载后执行 `sha256sum -c SHA256SUMS` 可检查完整性；COUL90 的独立校验值位于 `app/coulomb/`。输入或文档修改后，相关文件的哈希会相应变化。

## 实验比较与说明书

[说明书](docs/EM500核力与两核子态求解.pdf)第9章给出氘核比较表，以及 np、pp 九能量角分布。三个目录中的绘图程序直接读取已保存的数值数据。重新查询 EXFOR 或从论文提取数据属于独立步骤，所需来源和环境在各目录说明中列出。

物理求解器由 Fortran 实现。Python 用于部分绘图和文档排版；文档环境见 [PDF生成说明](docs/PDF生成说明.md)。修改正文后，在根目录运行：

```sh
python3 docs/build_pdf.py
```

## 物理范围与现有检验

- 计算对象为两核子束缚态与弹性散射。nd 散射、破裂和三体力属于其他项目中的三体计算。
- 氘核结合能参与过核力定参，复现它用于检查实现。D 态概率是依赖表示的诊断量。
- 氘核电荷半径采用说明书所列的有限核子尺寸修正；四极矩和磁矩采用单体电磁算符。完整二体电磁贡献属于进一步扩展内容。
- pp 采用当前点库仑、有效库仑强度和相同质子交换方案；完整精密电磁修正及显式介子产生通道处于当前模型范围之外。全角点库仑截面发散。
- 原 EM2003 构造报告的拟合数据低于290 MeV，300.2 MeV 比较用于检验拟合能区之外的描述。数值收敛与实验符合程度分别评估，现有偏差保留在说明书及比较资料中。

## 许可与引用

项目新增的接口、求解器、观测量和运行管理代码由 **sujun0106** 以 GNU GPLv3（GPL-3.0-only）发布，见 [LICENSE](LICENSE)。第三方组件和数据按各自来源处理，具体见 [第三方说明](THIRD_PARTY_NOTICES.md)。

EM500 文件取自 [Tic-tac 固定提交](https://github.com/seanbsm/Tic-tac/tree/30338ac0eaf968b6db2edfa8380410651a56267b)，保留源码署名和上游许可。使用本程序开展研究时，请引用 [Entem–Machleidt 2003](https://doi.org/10.1103/PhysRevC.68.041001)、对应数值算法以及使用的实验数据来源。
