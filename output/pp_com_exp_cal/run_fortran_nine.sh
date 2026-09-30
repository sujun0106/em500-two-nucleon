#!/usr/bin/env bash
# Reproducible nine-energy pp calculations, performed entirely by Fortran.
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
code=$(cd -- "$here/../.." && pwd)
make -C "$code" all
batch=${1:-"$here/fortran_runs/$(date -u +%Y%m%d_%H%M%S_UTC)"}
mkdir -p -- "$(dirname -- "$batch")"
mkdir -- "$batch"
batch=$(cd -- "$batch" && pwd)
sha256sum "$code/app/fewbody.x" > "$batch/executable.sha256"
date -u +%FT%TZ > "$batch/started_utc.txt"
for token in 3p037 9p69 14p16 50p06 98 142 239p9 270p8 286p7; do
  run="$batch/${token}MeV"
  mkdir -- "$run"
  cp -- "$here/fortran_inputs/pp_${token}MeV.nml" "$run/input.nml"
  echo "Starting pp ${token} MeV"
  (cd -- "$run" && "$code/app/fewbody.x" input.nml > run.log 2>&1)
  tail -n 3 "$run/run.log"
done
date -u +%FT%TZ > "$batch/completed_utc.txt"
printf 'fortran_runs/%s\n' "$(basename -- "$batch")" > "$here/latest_fortran_batch.txt"
echo "Completed batch: $batch"
