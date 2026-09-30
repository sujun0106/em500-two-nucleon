# Fortran EM500与九能量np实验角分布比较

- `np_exp_cal_angular_3x3.png` / `.pdf`：新的3×3实验—计算对比图。黑色实线为本次实际运行Fortran程序的结果；符号和误差棒为EXFOR实验数据。原纯实验图保留。
- 能量：3.1、10.04、25.8、50、96、137、150、200、300.2 MeV（入射中子实验室动能）。
- 每幅图横轴为质心系出射中子角度，纵轴为mb/sr，各子图纵轴刻度独立。
- 固定EM500两体核力，未重新拟合势参数，也未对理论曲线进行实验归一化。
- 动量Gauss网格96点，积分上限1000 MeV/c，总角动量J≤32；0—180°每0.25°输出一点，共721点。
- `relativistic`仅使用相对论运动学计算在壳动量，不意味着全相对论动力学。
- 300.2 MeV使用Franz 2000的111个实验点，按论文Table I的10.63 mb/sr绝对因子换算，额外3.1%共同归一化误差单列，未加进逐点误差棒。实验为20 MeV能量分箱，曲线为名义能量计算，未作能谱卷积。已移除284 MeV及其临时单位修正。

## 计算文件与复现

`fortran_inputs/`保存九个输入。`fortran_runs/`保存每次运行的输入快照、日志、角分布、S/K/T矩阵、分波、振幅和诊断；`latest_fortran_batch.txt`指向本图所用批次。

在本目录执行：

```sh
./run_fortran_nine.sh
python3 plot_np_comparison.py
```

第一步调用上两级code/app/fewbody.x（必要时先make），每次新建带UTC时间的结果批次，不覆盖原结果。计算全由Fortran完成；Python仅整理数据和绘图。直接重画只需第二步，不需要重复Fortran计算。服务器计算所需依赖仍以code/README.md为准。

`fortran_summary.csv`汇总总弹性截面与内部检查。`fortran_provenance.json`记录源码和可执行文件哈希、输入设置。`comparison_at_experimental_angles.csv`把0.25°计算网格线性插值到实验角度，提供逐点差值；此表不是在每个实验角度重新求解，也未进行束流展宽或探测器接受度卷积。原始计算采样点完整保留，图中的线直接连接这些采样点。

## 数值检查与解释边界

九组计算均正常结束。数值最大值如下：

- `angular_integral_relative_error`：1.2995e-15。
- `optical_relative_error`：3.24876e-15。
- `maximum_residual_relative`：3.44492e-14。
- `maximum_unitarity_error`：7.77156e-16。
- `J24_to_J32_max_relative_change_at_0_90_180_deg`：0.000464151。

内部残差、幺正性和角积分一致性不等于动量网格或分波已严格收敛。J24与J32的比较这里只覆盖0、90、180度，不能代替完整角区间的分波收敛研究；本轮未新增动量网格扫描。

图中实验覆盖范围有限；不对未测角度补实验点。这里比较名义单能计算与所选实验，不作能谱或角接受度卷积。未给出把所有误差视为独立的χ²结论。实NN-only模型未显式描述介子产生等非弹性通道，内部光学闭合不等于完整真实反应总截面。

## 上游原件

论文、论文页截图及NN-OnLine网页原件以[SOURCE_LINKS.md](SOURCE_LINKS.md#franz)中的原站链接提供；来源清单保留历史下载哈希。重画使用冻结数表与Fortran结果，无需原件。
