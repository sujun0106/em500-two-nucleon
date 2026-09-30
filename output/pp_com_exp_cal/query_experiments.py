#!/usr/bin/env python3
"""Optional EXFOR requery; frozen local data suffice for plotting."""
from pathlib import Path
import argparse,json,sys

def query_environment(parser, args):
    """Resolve optional external EXFOR inputs before changing saved data."""
    tool_root = args.workspace / 'tools/EXFOR' if args.workspace else None
    index = args.index or (tool_root / '.cache/exfor.sqlite' if tool_root else None)
    source = args.source_root or (tool_root / 'external/exfor-git/exforall' if tool_root else None)
    if index is None or source is None:
        parser.error('Requery is optional. Supply --index and --source-root (or --workspace). '
                     'Saved data and plot scripts work without an EXFOR database.')
    index, source = index.expanduser().resolve(), source.expanduser().resolve()
    if not index.is_file():
        parser.error(f'EXFOR SQLite index does not exist: {index}')
    if not source.is_dir():
        parser.error(f'EXFOR source directory does not exist: {source}')
    if args.tool_src:
        tool_src = args.tool_src.expanduser().resolve()
        if not (tool_src / 'exfor_tool/query.py').is_file():
            parser.error('--tool-src must contain exfor_tool/query.py')
        sys.path.insert(0, str(tool_src))
    return index, source


def arguments():
    parser = argparse.ArgumentParser(
        description='Optional EXFOR requery; offline plots use the saved query results.')
    parser.add_argument('--index', type=Path, help='External EXFOR SQLite index')
    parser.add_argument('--source-root', type=Path, help='External exforall source directory')
    parser.add_argument('--tool-src', type=Path,
                        help='Optional directory containing exfor_tool; default: bundled snapshot')
    parser.add_argument('--workspace', type=Path,
                        help='Compatibility shortcut for <workspace>/tools/EXFOR database paths')
    args = parser.parse_args()
    return parser, args

def main():
    root=Path(__file__).resolve().parent
    parser,args=arguments()
    index,source=query_environment(parser,args)
    from exfor_tool.query import execute_query
    from exfor_tool.render import write_table
    from exfor_tool.pp_angular_grid import extract_pp_cm
    config=json.loads((root/'panel_queries.json').read_text())
    for energy,path in zip(config['energies_MeV'],config['outcomes']):
        stem=f'{energy:g}MeV'.replace('.','p')
        request=json.loads((root/'queries'/f'request_{stem}.json').read_text())
        outcome=execute_query(request,index_path=index,source_root=source)
        if outcome['status']!='result':raise RuntimeError(outcome)
        (root/path).write_text(json.dumps(outcome,ensure_ascii=False,indent=2)+'\n')
        write_table(outcome,root/'queries'/f'raw_{stem}')
        (root/'queries'/f'prepared_panel_{stem}.json').write_text(json.dumps(extract_pp_cm(outcome),ensure_ascii=False,indent=2)+'\n')
        print(energy,len(outcome['records']))
if __name__=='__main__':main()
