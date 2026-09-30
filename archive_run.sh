#!/usr/bin/env bash
# Archive one completed run. The root input file stays editable for the next run.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
if (( $# )); then echo 'Usage: ./archive_run.sh' >&2; exit 2; fi
if [[ -d .run.lock && ${DEUTERON_AUTO_ARCHIVE:-0} != 1 ]]; then
  echo 'A run is active; wait for run.sh to finish before manual archiving.' >&2; exit 2
fi
status_file=run_status.txt
if [[ ! -f "$status_file" && -f deuteron_status.txt ]]; then status_file=deuteron_status.txt; fi
[[ -f "$status_file" ]] || { echo 'No completed run to archive.' >&2; exit 2; }
if ! mkdir .archive.lock 2>/dev/null; then echo 'Another archive operation is active.' >&2; exit 2; fi
trap 'rmdir .archive.lock 2>/dev/null || true' EXIT
value() { awk -F= -v key="$1" '{k=$1;gsub(/[[:space:]]/,"",k);if(k==key){v=$2;gsub(/[[:space:]]/,"",v);print v}}' "$status_file"; }
[[ $(value status) == completed ]] || { echo 'Run has not completed successfully; keeping all files in the root.' >&2; exit 2; }
# A wrapper interruption/failure must not be presented as a successful archive.
if [[ -f run_metadata.txt ]] && grep -Eq '^driver_status=(failed|interrupted)$' run_metadata.txt; then
  echo 'Driver recorded a failure/interruption; keeping the run for inspection.' >&2; exit 2
fi
finish_date=$(value finished_date)
finish_time=$(value finished_time)
finish_zone=$(value finished_zone)
[[ $finish_date =~ ^[0-9]{8}$ && $finish_time =~ ^[0-9]{6}(\.[0-9]+)?$ && $finish_zone =~ ^[+-][0-9]{4}$ ]] || {
  echo 'Missing or invalid completion timestamp.' >&2; exit 2;
}
iso="${finish_date:0:4}-${finish_date:4:2}-${finish_date:6:2} ${finish_time:0:2}:${finish_time:2:2}:${finish_time:4:2} $finish_zone"
stamp=$(date -u -d "$iso" +%Y%m%d_%H%M%S_UTC)
files=(run_input.nml "$status_file")
case "$(value calculation)" in
  deuteron_bound_state)
    effective_input=deuteron_input.nml
    files+=(deuteron_input.nml deuteron_summary.txt deuteron_wavefunction.csv)
    ;;
  np_scattering)
    effective_input=np_input.nml
    files+=(np_input.nml np_summary.txt np_partial_waves.csv np_on_shell.csv np_angular_distribution.csv np_j_convergence.csv)
    ;;
  pp_scattering)
    effective_input=pp_input.nml
    files+=(pp_input.nml pp_summary.txt pp_partial_waves.csv pp_nuclear_s.csv pp_angular_distribution.csv pp_j_convergence.csv)
    ;;
  *) echo 'Unknown calculation in status record; nothing archived.' >&2; exit 2 ;;
esac
for f in "${files[@]}"; do
  [[ -f "$f" && ! -L "$f" ]] || { echo "Missing regular run file: $f" >&2; exit 2; }
done
for f in deuteron_observables.dat deuteron_coordinate_wavefunction.csv deuteron_matrices.csv np_spin_amplitudes.csv pp_spin_amplitudes.csv run.log run_metadata.txt job.pid; do
  if [[ -e "$f" || -L "$f" ]]; then
    [[ -f "$f" && ! -L "$f" ]] || { echo "Not a regular run file: $f" >&2; exit 2; }
    files+=("$f")
  fi
done
mkdir -p output
suffix=0
while :; do
  if (( suffix==0 )); then destination="output/$stamp"; else printf -v destination 'output/%s_%02d' "$stamp" "$suffix"; fi
  if mkdir -- "$destination" 2>/dev/null; then break; fi
  [[ -e "$destination" ]] || { echo "Cannot create $destination" >&2; exit 2; }
  suffix=$((suffix+1))
done
# Copy, compare, and only then remove the known runtime files. If copying fails,
# root originals remain available. Never move an arbitrary glob of root files.
archived=()
for f in "${files[@]}"; do
  target=$f
  if [[ $f == run_input.nml ]]; then target=input.nml; fi
  cp -p -- "$f" "$destination/$target"
  cmp -s -- "$f" "$destination/$target"
  archived+=("$target")
done
{
  printf 'archive_status=completed\n'
  printf 'calculation_finished_local=%s %s %s\n' "$finish_date" "$finish_time" "$finish_zone"
  printf 'directory_time_basis=calculation completion converted to UTC\n'
  printf 'archived_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'input.nml is the original snapshot actually read by the calculation.\n'
  printf '%s contains the effective interpreted parameters.\n' "$effective_input"
} > "$destination/archive_info.txt"
archived+=(archive_info.txt)
(cd -- "$destination" && sha256sum -- "${archived[@]}" > SHA256SUMS && sha256sum -c SHA256SUMS >/dev/null)
rm -- "${files[@]}"
printf 'Archived: %s\n' "$destination"
