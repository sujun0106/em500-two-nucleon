# 第三方来源与许可说明

项目新增的接口、两体求解器、观测量和运行脚本由 sujun0106 以 GNU GPLv3（GPL-3.0-only）发布，完整许可见根目录 LICENSE。下表分别说明第三方组件的来源与分发方式。

| 组成 | 来源记录 | 公开版处理 |
|---|---|---|
| EM500 Fortran核力实现 | app/provenance/upstream.json、源码署名及app/em500/LICENSE | 保留Tic-tac固定提交30338ac0eaf968b6db2edfa8380410651a56267b中的实现、原署名与上游GPLv3文本 |
| LAPACK/BLAS子集 | app/lapack/UPSTREAM.json、LICENSE、SHA256SUMS | 按上游BSD-3-Clause条件保留版权、许可及免责声明 |
| COUL90 | app/coulomb/UPSTREAM.json、README.md | 仓库保存来源和固定校验值；使用者通过fetch_external.sh从原站获取，文件在本地被Git忽略 |
| 实验数值数据 | output三个比较目录中的来源清单 | 保留数值、单位、作者、数据库和论文出处；数据权利按相应来源处理 |
| 第三方论文全文、论文页截图与原网页全文 | 各比较目录的SOURCE_LINKS.md及历史下载清单 | 公开包提供原始来源链接，研究归档中的全文副本独立保留 |
| 本项目制作的比较图、计算输入和结果 | output中的Fortran运行记录和绘图脚本 | 随公开包提供以便复现，引用实验结果时应同时引用原数据来源 |

COUL90作者为A. R. Barnett，原站为 https://www.fresco.org.uk/functions/barnett/index.htm 。原文件和下载页面的再分发许可尚待明确，因此本公开包通过来源链接和校验脚本准备本地依赖。项目GPLv3声明的适用范围为项目代码；COUL90的使用条件仍以原来源为准。

v0.5.2调整分发方式及文档。固定核力配置、Fortran求解算法和已保存的物理数值均沿用v0.5.1。
