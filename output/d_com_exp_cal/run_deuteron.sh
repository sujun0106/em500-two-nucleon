#!/usr/bin/env bash
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
code=$(cd -- "$here/../.." && pwd)
make -C "$code" all
run="$here/fortran_runs/$(date -u +%Y%m%d_%H%M%S_UTC)"
mkdir -p -- "$(dirname -- "$run")"
mkdir -- "$run"
cp -- "$here/input.nml" "$run/input.nml"
sha256sum "$code/app/fewbody.x" > "$run/executable.sha256"
cp -- "$code/VERSION" "$run/VERSION"
(cd -- "$code" && find app -type f \( -iname '*.f90' -o -iname '*.f' -o -name '*.mk' \) -print0 | sort -z | xargs -0 sha256sum) > "$run/source_files.sha256"
(cd -- "$code" && sha256sum Makefile) >> "$run/source_files.sha256"
(cd -- "$run" && "$code/app/fewbody.x" input.nml > run.log 2>&1)
printf 'fortran_runs/%s\n' "$(basename -- "$run")" > "$here/latest_run.txt"
cat "$run/deuteron_observables.dat"
