# 重新生成《EM500核力与两核子态求解》PDF

本文对应公开代码包0.5.2、主说明书v0.3，更新日期为2026-09-30。`build_pdf.py` 直接读取正文与已有自制图，完成Markdown转换和XeLaTeX排版；构建输入全部采用包内相对路径。物理求解、实验数据查询及 COUL90获取是独立的工作步骤。仅生成说明书时，准备下列排版环境即可。

## 输入与输出

- `docs/EM500核力与两核子态求解.md`：正文主源；默认先转换为同名 `.tex`。
- `docs/说明书素材/example_outputs.pdf`（也可仅保留 `.png`）：第7章三个运行范例的已保存图。
- `output/np_com_exp_cal/np_exp_cal_angular_3x3.pdf`（或 `.png`）：第9章np九宫格，末组为300.2 MeV。
- `output/pp_com_exp_cal/pp_exp_cal_angular_3x3.pdf`（或 `.png`）：第9章pp九宫格，对数纵轴。
- `output/d_com_exp_cal/`：氘核实验比较与来源记录；当前手册表格直接保存在正文中，排版时不重新生成它。

图片优先嵌入同名矢量PDF，不存在时读取Markdown中指定的PNG。最终PDF和TeX保存在 `docs/`。已有 `docs/EM500核力与两核子态求解.tex` 也可直接编译；这一路径不需要Markdown解析包。正文的历史核验记录指向当时的运行证据，不构成排版依赖。

## 环境

Python 3、XeLaTeX、CTeX与Fandol中文字体、DejaVu字体。TeX还需 `geometry`、`amsmath`、`amssymb`、`bm`、`graphicx`、`longtable`、`booktabs`、`array`、`enumitem`、`fvextra`、`xcolor`、`xurl`、`fancyhdr`、`needspace`、`hyperref` 等常见宏包。Debian/Ubuntu上可安装：

```sh
sudo apt-get install texlive-xetex texlive-lang-chinese texlive-latex-extra fonts-dejavu-core
python3 -m pip install -r docs/requirements-pdf.txt
```

文档排版使用上述环境；Fortran计算的依赖由主手册第7.1节单独说明。Markdown依赖也可在Python虚拟环境中安装。

## 从Markdown更新PDF

在代码包根目录执行：

```sh
python3 docs/build_pdf.py
```

也可从任意目录给出脚本路径运行。脚本按自身位置找文件，默认在系统临时目录进行三次XeLaTeX编译，成功后更新PDF并清理临时文件。

需要保存日志时：

```sh
python3 docs/build_pdf.py --build-dir /tmp/em500-manual-build
```

## 仅用现有TeX更新PDF

```sh
python3 docs/build_pdf.py --from-tex
```

或不使用Python，直接在 `docs/` 编译：

```sh
cd docs
xelatex -interaction=nonstopmode -halt-on-error EM500核力与两核子态求解.tex
xelatex -interaction=nonstopmode -halt-on-error EM500核力与两核子态求解.tex
xelatex -interaction=nonstopmode -halt-on-error EM500核力与两核子态求解.tex
```

直接调用XeLaTeX会在 `docs/` 留下 `.aux`、`.log`、`.toc`、`.out` 等中间文件。推荐使用脚本管理这些文件。

主手册的九宫格直接嵌入包内已有实验与计算比较图，保留原数值、误差棒及图注。公开版通过各实验比较目录的 `SOURCE_LINKS.md` 提供第三方论文和网页入口；这些第三方全文、截图及网页快照均属于排版输入之外的资料。COUL90获取脚本也属于物理构建流程，PDF生成使用独立的文档依赖。
