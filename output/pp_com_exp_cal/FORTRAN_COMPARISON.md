# pp Fortran计算与验证说明

全部九能量已完成实际Fortran计算，模型是固定EM500的pp核力加点库仑相互作用，包含两个相同质子的振幅反对称化。未调整核力或曲线归一化来贴合实验。未另加两光子、真空极化、磁矩等电磁修正及显式非弹性通道。

使用160点动量网格，pmax=1000 MeV/c，J≤32，库仑匹配半径16 fm，径向积分192点；每能量输出0.25—179.75度719点。所有细节在各运行的input.nml及pp_summary.txt中。

## 文件与复现

- `run_fortran_nine.sh`：从现有code编译并运行九能量，自动新建批次；物理计算不依赖Python。
- `fortran_inputs/`：九份输入；`fortran_runs/`：输入快照、日志、角分布、核S矩阵、振幅及诊断。
- `latest_fortran_batch.txt`：当前比较图所用批次。
- `fortran_summary.csv`、`calculation_validation.json`：诊断与验证。
- `comparison_at_experimental_angles.csv`：将0.25度计算网格线性插值到实验角度的差值表。未对实验做插值，也未做束流/角接受度卷积。
- `numerical_sensitivity.csv`：低/高能完整角分布数值参数变化检查。
- `independent_PWA93_comparison.csv`：独立经验相移分析的45点参考，仅作量级/形状检查，未画入实验图，不属于同EM500核力的作者结果。

```sh
./run_fortran_nine.sh
python3 plot_pp_comparison.py
```

若只重画，运行第二行。独立实验图由`plot_pp_3x3.py`重画。Python只处理数据及作图；服务器上的Fortran核力、求解、振幅计算仍使用同一个`app/fewbody.x`。

## 验证与边界

原氘核/np回归通过；新增pp原核力电荷接口、自由Bessel/Coulomb函数、α→0桥接、纯库仑Mott截面、库仑—核干涉相位、Pauli张量交换及错误输入检查通过。既有核力原始两文件没有改动。

九个计算最大诊断：

- `maximum_residual_relative`：4.03156e-15。
- `maximum_reciprocity_error_relative`：1.61563e-15。
- `maximum_unitarity_error`：9.99201e-16。
- `theta_reflection_relative_error`：1.25434e-15。
- `J24_to_J32_max_relative_change_at_10_30_60_90_150deg`：4.53988e-05。

完整719角度敏感性：

- `high_N112`，286.7 MeV：最大相对变化2.01386e-14。
- `high_R12`，286.7 MeV：最大相对变化3.14648e-06。
- `low_R12`，3.037 MeV：最大相对变化2.01684e-06。

上述检查支持本轮画图所需的数值稳定性，不等于所有参数和模型误差都为零。J24→32的诊断只在所列五个角度统计。纯库仑在有限pmax=1000下仍有最大约1.21×10⁻⁵的|S_C−I|误差；提高动量上限会降低它，不能仅靠增加动量点数消除。它与幺正性误差不同。

实验误差定义见DATA_SOURCES.md。9.69 MeV记录给出概然误差；多组有共同归一化误差，因此未把全部误差视作独立1σ进行χ²判决。较高能区计算与实验仍有偏差。当前验证不能把这种差异唯一归因于实验或EM500强核力；完整电磁处理、模型实现的外部同核力验证及实验接受度处理仍有进一步核查空间。

库仑在0/180度发散，本图纵轴按实测角域设置，端点尖峰可能超出显示范围。输出的是微分截面，未宣称有限的全角库仑总截面。

数值方法见[《EM500核力与两核子态求解》](../../docs/EM500核力与两核子态求解.pdf)第5章，数值检查见第8章。源码、COUL90来源和校验文件随code迁移。项目历史完整审计存于原项目的 `runs/pp_extension_2026-09-29/`，不是此发布包运行或重画的依赖。

纵轴已按用户要求改为对数坐标（截面数值及单位mb/sr不变）；各子图范围独立，PNG/PDF与重画代码同步更新。
