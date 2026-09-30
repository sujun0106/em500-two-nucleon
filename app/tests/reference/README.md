# 氘核迁移参考

本目录只含文本，不含 Python 程序。`make check` 可以用纯 Fortran 读取数值表。

- `reference_bound_states.dat`：前两行为注释，第三行为数据条数；随后为 `N pmax_MeV energy_MeV probability_S probability_D kinetic_MeV potential_MeV`。八组数值直接提取自历史氘核运行，不是新 Fortran 程序自产的“标准答案”。
- `momentum_wavefunction_N100_pmax1000.csv`：逐字节复制历史 N=100、pmax=1000 MeV/c 波函数。列为 `p, weight, psi_S, psi_D, p² psi_S², p² psi_D²`。波函数归一化为 `sum_i weight_i p_i² (psi_S(i)² + psi_D(i)²) = 1`。
- `historical_convergence.csv`：完整保留历史收敛表，包括当时的残差、厄米性指标和运行时间。运行时间不作迁移验收标准。
- `reference_potential.dat`：前三行为注释，第四行为数据条数；随后为 `pout pin V_SS V_SD V_DS V_DD`。由原 Python API 调用原已编译 EM500 Fortran 库生成。含不对称动量及交换方向，避免只检验对角矩阵元而遗漏 S/D 行列互换。
- `reference_manifest.json`：历史代码版本、源文件路径、SHA256、势矩阵约定和限制。
- `SHA256SUMS`：上述文件校验值，可执行 `cd reference && sha256sum -c SHA256SUMS`。

势矩阵元单位为 MeV^-2，动量为 MeV/c，积分测度为 `dp p²`，不另乘 `2/pi`。矩阵第一下标是出射轨道角动量，第二下标是入射轨道角动量。原核力输出的 `pot(4), pot(6), pot(5), pot(3)` 分别对应 `SS, SD, DS, DD`。采用原 `radial_rephased` 相位约定；比较整个波函数时允许同一个整体正负号，不能只翻转 D 波符号。

数值积分上限 `pmax=1000 MeV/c` 不是调节器尺度 `500 MeV`。本上游实现在任意一个动量大于 1000 MeV/c 时返回零；边界参考用于检查保留该实现行为。

这些数据检验语言迁移的一致性，不是独立核力物理验证。结合能曾参与相互作用标定；D 波概率依赖表示，不是直接实验可观测量。

历史计算 Tool 版本：`57713f32c18f0bebb840c027e081c096fdac22cc`（`v0.1.0-alpha.1`）。势样本生成时当前工作树版本：`026aef24f4cc0d57fbabf6f4f9f5c4e23fda58a5`；所用旧库哈希与历史记录一致。参考生成脚本留在项目运行记录 `runs/fortran_deuteron_2026-09-29/reference_audit/prepare_reference.py`，服务器运行不需要它。
