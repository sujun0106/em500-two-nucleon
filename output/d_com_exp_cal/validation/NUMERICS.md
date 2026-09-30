# 氘核坐标波函数与观测量数值核验

本目录是开发审计，不是Fortran发布包的运行依赖。正式实现为 `code/app/deuteron_observables.f90`；根目录程序由主任务整合。这里未修改EM500核力参数或结合能拟合。

## 定义和相位

动量波函数归一化为 ∑L∫dp p²ψL(p)²=1，p以MeV计。坐标波函数：

`u_L(r)=sqrt(2/pi) r/(hbarc)^(3/2) ∫dp p² j_L(pr/hbarc) ψ_L(p)`，L=0、2。

现有EM500 `radial_rephased` 约定满足 `psi_plain,L=i^L psi_rad,L`，Fourier核展开再提供i^L，因此L=0、2均为正号Hankel变换；不能单独给D分量加负号。独立解析Gaussian S+D测试检验了相对符号及四极矩干涉项；非通过“匹配实验四极矩”挑选符号。

r是np相对距离，等质量质心近似R=r/2。每核子球平均密度为 `rho1(R)=[u(2R)^2+w(2R)^2]/(2*pi*R^2)`，∫4πR²rho1dR=1；总核子密度为2rho1，归一为2。R=0用S波解析斜率求极限。它是点核子密度，不是实验电荷密度，未进行质子/中子电荷形状折叠。

`r_point²=(1/4)∫dr r²(u²+w²)`；`Q_IA=(1/20)∫dr r² w(sqrt(8)u-w)`。代码采用原模型hbarc=197.327 MeV fm。有限坐标网格的波函数不再次归一化，以保留截尾/变换误差诊断。

A_S和eta取固定gamma=sqrt(mB)/hbarc的20–25fm线性最小二乘拟合，`u=A_S exp(-gamma r)`，`w=A_D exp(-gamma r)(1+3/(gamma r)+3/(gamma r)^2)`，eta=A_D/A_S。

## 核验结果

- 独立解析S与S+D Gaussian变换：逐点u,w、Parseval归一、半径、原点密度Jacobian、四极矩均通过1e-11量级容差。测试源码在 `code/app/tests/check_deuteron_observables.f90`。
- 扫描N_p=100、160、192、240，pmax=1000MeV；rmax=30、40、50、60fm，dr=0.01fm，记录于 `convergence.csv`。
- 最终N_p=192、rmax=60fm、6000坐标区间：norm=.9999999999985604，r_point=1.974496510325755fm，Q_IA=.274951137566983fm²。
- N_p=160→240，半径变化<3e-13fm；rmax=50→60，半径变化2.51e-8fm，Q变化1.80e-9fm²。
- dr=.02/.01/.005fm比较，半径/Q变化约1e-15，见 `radial_spacing_differences.json`。
- 独立动量空间公式 `<r_rel²>=hbarc² sum_L ∫dp[p² (dψ_L/dp)²+L(L+1)ψ_L²]`，采用Gauss-Legendre谱微分得r_point=1.974496510677491fm，与坐标结果差3.52e-10fm，见 `momentum_radius_audit.py/json`。
- 20–25fm拟合A_S=.884302876502fm^(-1/2)，eta=.025580730316。18–23、22–27、25–30fm窗口的eta分别.025580784、.025580746、.025580949；约2.2e-7的拟合窗口变化应与物理模型误差区分。窗口内点系数仍有小幅变化，未将其伪装成精确常数，全部范围在 `tail_windows.json`。

## 物理精度界限

数值收敛不能保证电磁可观测量的完备性。Q和磁矩是单体脉冲近似；电荷半径估计只包含核子尺寸及Darwin–Foldy项，相对论及二体电磁修正未完整实现。波函数和D态概率依赖表示，不作为直接实验观测量。渐近参数是从数据分析提取的量，也不等同于直接测量每个r处的波函数。

## 复现

在 `code/` 根目录执行 `make check` 运行正式保留的验收测试；解析变换测试源码为 `app/tests/check_deuteron_observables.f90`。这里保存的收敛表是历史扫描结果，不依赖未随包提供的开发脚本。服务器正式运行仍使用根目录 Makefile 及内置 LAPACK 模式。momentum_radius_audit.py为独立验证脚本，依赖NumPy；正式观测量计算不依赖Python。
