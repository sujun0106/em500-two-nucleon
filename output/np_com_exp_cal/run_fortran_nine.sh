#!/usr/bin/env bash
# Run the unmodified production Fortran solver in separate result directories.
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
code=$(cd -- "$here/../.." && pwd)
make -C "$code" all
batch=${1:-"$here/fortran_runs/$(date -u +%Y%m%d_%H%M%S_UTC)"}
mkdir -p -- "$(dirname -- "$batch")"
mkdir -- "$batch"  # Refuse overwriting an existing batch.
batch=$(cd -- "$batch" && pwd)
sha256sum "$code/app/fewbody.x" > "$batch/executable.sha256"
date -u +%FT%TZ > "$batch/started_utc.txt"
for token in 3p1 10p04 25p8 50 96 137 150 200 300p2; do
  run="$batch/${token}MeV"
  mkdir -- "$run"
  cp -- "$here/fortran_inputs/np_${token}MeV.nml" "$run/input.nml"
  echo "Starting ${token} MeV"
  (cd -- "$run" && "$code/app/fewbody.x" input.nml > run.log 2>&1)
  tail -n 4 "$run/run.log"
done
date -u +%FT%TZ > "$batch/completed_utc.txt"
# Keep the pointer relocatable, including an explicitly selected batch folder.
common=$here
upward=
while [[ "$batch" != "$common/"* ]]; do
  common=${common%/*}
  upward="../$upward"
done
printf '%s%s\n' "$upward" "${batch#"$common/"}" > "$here/latest_fortran_batch.txt"
echo "Completed batch: $batch"
