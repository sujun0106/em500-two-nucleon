#!/usr/bin/env bash
# Generic Linux background submission. This is not a Slurm/PBS queue script.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
[[ -x ./app/fewbody.x ]] || { echo 'Compile first: make' >&2; exit 2; }
# Reject invalid arguments before detaching, so the user sees the reason.
check_args=("$@")
if [[ ${check_args[0]:-} == --no-archive ]]; then check_args=("${check_args[@]:1}"); fi
if (( ${#check_args[@]} > 1 )); then echo 'Usage: ./submit.sh [--no-archive] [input.nml]' >&2; exit 2; fi
[[ -f ${check_args[0]:-input.nml} ]] || { echo "Input not found: ${check_args[0]:-input.nml}" >&2; exit 2; }
if [[ -d .run.lock ]]; then echo 'A run is already active, or its lock needs inspection.' >&2; exit 2; fi
source ./app/runtime_files.sh
for f in "${runtime_files[@]}"; do
  if [[ -e "$f" || -L "$f" ]]; then echo "Existing runtime file: $f; archive or preserve it before submission." >&2; exit 2; fi
done
nohup env DEUTERON_BACKGROUND=1 ./run.sh "$@" >/dev/null 2>&1 < /dev/null &
job_pid=$!
printf 'Background PID: %s\n' "$job_pid"
echo 'Runtime log: run.log; successful runs are archived under output/<completion-time>/.'
