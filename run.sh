#!/usr/bin/env bash
# Run in the package root, then archive using the program's completion time.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
auto_archive=1
if [[ ${1:-} == --no-archive ]]; then auto_archive=0; shift; fi
if (( $# > 1 )); then echo 'Usage: ./run.sh [--no-archive] [input.nml]' >&2; exit 2; fi
input_file=${1:-input.nml}
[[ -f "$input_file" ]] || { echo "Input not found: $input_file" >&2; exit 2; }
[[ -x ./app/fewbody.x ]] || { echo 'Compile first: make' >&2; exit 2; }
# Only one run uses the fixed root filenames. No old result is overwritten.
if ! mkdir .run.lock 2>/dev/null; then
  echo 'A run lock exists. Check the active job before removing a stale .run.lock.' >&2
  exit 2
fi
trap 'rmdir .run.lock 2>/dev/null || true' EXIT
source ./app/runtime_files.sh
for f in "${runtime_files[@]}"; do
  if [[ -e "$f" || -L "$f" ]]; then
    echo "Existing runtime file: $f. Archive a completed run, or move failed-run files aside first." >&2
    exit 2
  fi
done
if [[ ${DEUTERON_BACKGROUND:-0} == 1 ]]; then printf '%s\n' "$$" > job.pid; fi
{
  printf 'started_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'input_source=%s\n' "$input_file"
  printf 'driver_pid=%s\n' "$$"
  printf 'package_version=%s\n' "$(cat VERSION)"
  printf 'executable_sha256_at_start=%s\n' "$(sha256sum app/fewbody.x | awk '{print $1}')"
} > run_metadata.txt
child=''
stop_run() {
  local signal=$1 code=$2
  trap - INT TERM HUP
  if [[ -n "$child" ]]; then
    kill -TERM "$child" 2>/dev/null || true
    wait "$child" 2>/dev/null || true
  fi
  printf 'driver_status=interrupted\nsignal=%s\n' "$signal" >> run_metadata.txt
  echo 'Interrupted; runtime files remain in the root directory.' >&2
  exit "$code"
}
trap 'stop_run INT 130' INT
trap 'stop_run TERM 143' TERM
trap 'stop_run HUP 129' HUP
./app/fewbody.x "$input_file" > run.log 2>&1 &
child=$!
if wait "$child"; then
  child=''
  printf 'driver_status=completed\n' >> run_metadata.txt
  cat run.log
  if (( auto_archive )); then DEUTERON_AUTO_ARCHIVE=1 ./archive_run.sh; else echo 'Completed. Run ./archive_run.sh when ready.'; fi
else
  rc=$?
  child=''
  printf 'driver_status=failed\nexit_code=%s\n' "$rc" >> run_metadata.txt
  cat run.log >&2
  echo 'Calculation failed; runtime files remain in the root directory.' >&2
  exit "$rc"
fi
