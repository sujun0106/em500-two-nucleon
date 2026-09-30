#!/usr/bin/env python3
"""Build the bundled EM500 manual without running physics or reading project runs/.

Dependencies: Python 3, markdown-it-py (unless --from-tex), XeLaTeX and fonts.
Run from any working directory; paths are relative to this script.
"""
from pathlib import Path
from urllib.parse import unquote
import argparse
import re
import shutil
import subprocess
import tempfile


def render_markdown(docs):
    try:
        from markdown_it import MarkdownIt
    except ImportError as exc:
        raise SystemExit('Install Markdown support: python3 -m pip install -r docs/requirements-pdf.txt; or use --from-tex.') from exc
    source=docs/'EM500核力与两核子态求解.md'
    text=source.read_text()
    maths=[]
    def protect(m,display):
     key=f'ZXMATH{len(maths):05d}XZ'
     maths.append((m.group(1),display))
     return ('\n\n'+key+'\n\n') if display else key
    # protect TeX before CommonMark interprets backslashes and pipes
    text=re.sub(r'\\\[(.*?)\\\]',lambda m:protect(m,True),text,flags=re.S)
    text=re.sub(r'\\\((.*?)\\\)',lambda m:protect(m,False),text,flags=re.S)
    md=MarkdownIt('commonmark').enable('table')
    tokens=md.parse(text)
    trans={'\\':r'\textbackslash{}','{':r'\{','}':r'\}','$':r'\$','&':r'\&','#':r'\#','%':r'\%','_':r'\_','^':r'\textasciicircum{}','~':r'\textasciitilde{}'}
    def escape(s):return ''.join(trans.get(c,c) for c in s)
    def rawmath(idx):
     tex,display=maths[int(idx)]
     if not display:return '\\('+tex+'\\)'
     # Keep numbers outside optional size reduction. Equations use the same font
     # normally; only over-wide equations are scaled to 92% of the text width.
     tags=re.findall(r'\\tag\{([^}]+)\}',tex)
     tex=re.sub(r'\\tag\{[^}]+\}','',tex).strip()
     tag=('\\tag{'+tags[0]+'}') if tags else ''
     return '\n\\begin{equation*}\\FitEquation{'+tex+'}'+tag+'\\end{equation*}\n'
    def text_render(s):
     chunks=re.split(r'ZXMATH(\d{5})XZ',s)
     return ''.join(rawmath(c) if i%2 else escape(c) for i,c in enumerate(chunks))
    def code(s):
     out=escape(s)
     # Permit wrapping at separators inside paths and output/parameter names.
     out=out.replace(r'\_',r'\_\allowbreak{}').replace('/',r'/\allowbreak{}').replace('.',r'.\allowbreak{}').replace('=',r'=\allowbreak{}')
     return '{\\ttfamily\\small '+out+'}'
    def inline(children):
     out=[]
     for t in children or []:
      if t.type=='text':out.append(text_render(t.content))
      elif t.type=='code_inline':out.append(code(t.content))
      elif t.type=='strong_open':out.append(r'\textbf{')
      elif t.type=='em_open':out.append(r'\emph{')
      elif t.type in ('strong_close','em_close'):out.append('}')
      elif t.type=='softbreak':out.append('\n')
      elif t.type=='hardbreak':out.append(r'\\'+'\n')
      elif t.type=='link_open':out.append('\\href{'+escape(t.attrGet('href'))+'}{')
      elif t.type=='link_close':out.append('}')
      elif t.type=='image':
       asset=Path(unquote(t.attrGet('src')))
       if asset.is_absolute():raise RuntimeError('Image paths must be relative: '+str(asset))
       preferred=asset.with_suffix('.pdf')
       if (docs/preferred).is_file():asset=preferred
       if not (docs/asset).is_file():raise FileNotFoundError(docs/asset)
       out.append('\n\\begin{center}\\includegraphics[width=\\linewidth,height=.72\\textheight,keepaspectratio]{'+asset.as_posix()+'}\\end{center}\n')
      else:raise RuntimeError(('inline',t.type))
     return ''.join(out)
    output=[];i=0;in_item=False
    while i<len(tokens):
     t=tokens[i]
     if t.type=='heading_open':
      level=int(t.tag[1]);title=tokens[i+1].content
      if level==1:
       output.append('\\begin{titlepage}\\centering\\vspace*{28mm}\n{\\Huge\\bfseries '+escape(title)+'\\par}\\vspace{12mm}\n{\\Large 程序说明书与操作教程\\par}\\vspace{8mm}\n{\\large 文档 v0.3 · 代码 v0.5.2\\par}\\vspace{8mm}\n2026年9月30日\\vfill\n\\begin{minipage}{.86\\textwidth}\\large\n从物理模型、连续方程与离散求解，到源码、输入和结果。\\par\\medskip\n包含氘核、np散射、pp散射三个逐步操作范例。\\end{minipage}\\vspace*{20mm}\\end{titlepage}\n\\tableofcontents\\clearpage\n')
      elif level==2:
       output.append('\\clearpage\\section*{'+inline(tokens[i+1].children)+'}\\addcontentsline{toc}{section}{'+escape(title)+'}\n')
      elif level==3:
       output.append('\\Needspace{6\\baselineskip}\\subsection*{'+inline(tokens[i+1].children)+'}\\addcontentsline{toc}{subsection}{'+escape(title)+'}\n')
      else:output.append('\\subsubsection*{'+inline(tokens[i+1].children)+'}\n')
      i+=3;continue
     if t.type=='table_open':
      rows=[];row=[];i+=1
      while tokens[i].type!='table_close':
       x=tokens[i]
       if x.type=='tr_open':row=[]
       elif x.type=='inline':row.append(inline(x.children))
       elif x.type=='tr_close':rows.append(row)
       i+=1
      n=len(rows[0]); widths=([.29,.34,.37] if n==3 else [.25]*4 if n==4 else [.38,.62] if n==2 else [1/n]*n)
      spec=''.join('>{\\raggedright\\arraybackslash}p{'+str(round((.94 if n>=4 else .96 if n==3 else .975)*w,4))+'\\linewidth}' for w in widths)
      output.append('\n{\\small\\setlength{\\tabcolsep}{3pt}\\renewcommand{\\arraystretch}{1.28}\n\\begin{longtable}{@{}'+spec+'@{}}\\toprule\n')
      output.append(' & '.join('\\textbf{'+c+'}' for c in rows[0])+r'\\\midrule'+'\n\\endfirsthead\n')
      output.append(' & '.join('\\textbf{'+c+'}' for c in rows[0])+r'\\\midrule'+'\n\\endhead\n')
      for row in rows[1:]:output.append(' & '.join(row)+r'\\'+'\n')
      output.append('\\bottomrule\\end{longtable}}\n')
     elif t.type=='paragraph_open':pass
     elif t.type=='paragraph_close':output.append('\n\n')
     elif t.type=='inline':output.append(inline(t.children))
     elif t.type=='fence':output.append('\n\\begin{Verbatim}[samepage=true,breaklines=true,breakanywhere=true,fontsize=\\small,frame=leftline,framesep=3mm,rulecolor=\\color{bluegray}]\n'+t.content+'\\end{Verbatim}\n')
     elif t.type=='bullet_list_open':output.append('\\begin{itemize}\n')
     elif t.type=='bullet_list_close':output.append('\\end{itemize}\n')
     elif t.type=='ordered_list_open':
      start=t.attrGet('start') or 1
      output.append('\\begin{enumerate}[start='+str(start)+']\n')
     elif t.type=='ordered_list_close':output.append('\\end{enumerate}\n')
     elif t.type=='list_item_open':output.append('\\item ')
     elif t.type=='list_item_close':output.append('\n')
     elif t.type=='html_block' and t.content.strip()=='<!-- pdf-pagebreak -->':output.append('\n\\clearpage\n')
     elif t.type=='hr':output.append('\\medskip\\hrule\\medskip\n')
     else:raise RuntimeError(('block',t.type))
     i+=1
    preamble=r'''\documentclass[UTF8,a4paper,zihao=5,fontset=fandol]{ctexart}
    \usepackage[margin=21mm,headheight=15pt]{geometry}
    \usepackage{amsmath,amssymb,bm,graphicx,longtable,booktabs,array,enumitem,fvextra,xcolor,xurl,fancyhdr,needspace}
    \usepackage{hyperref}
    \setmainfont{DejaVu Serif}
    \setsansfont{DejaVu Sans}
    \setmonofont[Scale=.88]{DejaVu Sans Mono}
    \definecolor{bluegray}{HTML}{285B7A}
    \hypersetup{colorlinks=true,linkcolor=bluegray,urlcolor=bluegray,pdftitle={EM500核力与两核子态求解},pdfauthor={含时量子多体调研项目}}
    \setlength{\parindent}{0pt}
    \setlength{\parskip}{.5em}
    \setlength{\emergencystretch}{3em}
    \setlist{nosep,leftmargin=2em}
    \setcounter{tocdepth}{2}
    \renewcommand{\baselinestretch}{1.14}
    \pagestyle{fancy}\fancyhf{}
    \fancyhead[L]{\small EM500核力与两核子态求解}
    \fancyhead[R]{\small v0.3 / code v0.5.2}
    \fancyfoot[C]{\thepage}
    \newsavebox{\EquationBox}
    \newcommand{\FitEquation}[1]{\sbox{\EquationBox}{$\displaystyle #1$}%
    \ifdim\wd\EquationBox>.92\linewidth\resizebox{.92\linewidth}{!}{\usebox{\EquationBox}}\else\usebox{\EquationBox}\fi}
    \allowdisplaybreaks
    \begin{document}
    '''
    tex=docs/'EM500核力与两核子态求解.tex'
    tex.write_text(preamble+''.join(output)+'\n\\end{document}\n')
    return tex, len(maths)


def compile_pdf(docs, tex, build):
    executable = shutil.which('xelatex')
    if executable is None:
        raise SystemExit('xelatex was not found. See docs/PDF生成说明.md for required TeX packages and fonts.')
    build.mkdir(parents=True, exist_ok=True)
    for run in range(3):
        completed = subprocess.run(
            [executable, '-interaction=nonstopmode', '-halt-on-error',
             f'-output-directory={build}', tex.name],
            cwd=docs, capture_output=True, text=True,
        )
        (build / f'pass{run + 1}.txt').write_text(completed.stdout + completed.stderr)
        if completed.returncode:
            print((completed.stdout + completed.stderr)[-6000:])
            raise RuntimeError('XeLaTeX failed. Re-run with --build-dir PATH to retain full logs.')
    built = build / (tex.stem + '.pdf')
    target = docs / built.name
    if built.resolve() != target.resolve():
        shutil.copy2(built, target)
    log = (build / (tex.stem + '.log')).read_text(errors='replace')
    markers = ('Missing character:', 'Overfull', 'Undefined', 'LaTeX Warning:')
    issues = [line for line in log.splitlines() if any(marker in line for marker in markers)]
    (build / 'pdf_issues.txt').write_text('\n'.join(issues) + ('\n' if issues else ''))
    print(f'PDF built: {target} ({target.stat().st_size} bytes); layout warnings: {len(issues)}')
    if issues:
        print('\n'.join(issues[:20]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--from-tex', action='store_true', help='Compile the existing .tex without Markdown/Python package dependencies.')
    parser.add_argument('--build-dir', type=Path, help='Retain TeX intermediate files and logs here; default: a temporary directory.')
    args = parser.parse_args()
    docs = Path(__file__).resolve().parent
    if args.from_tex:
        tex = docs / 'EM500核力与两核子态求解.tex'
        if not tex.is_file():
            raise SystemExit(f'Missing TeX source: {tex}')
    else:
        tex, count = render_markdown(docs)
        print(f'Markdown converted: {count} math expressions.')
    if args.build_dir is not None:
        compile_pdf(docs, tex, args.build_dir.expanduser().resolve())
    else:
        with tempfile.TemporaryDirectory(prefix='em500-manual-') as directory:
            compile_pdf(docs, tex, Path(directory))


if __name__ == '__main__':
    main()
