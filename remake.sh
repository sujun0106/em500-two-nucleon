#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
# Clean only generated compilation products directly inside app/.
if (( $# )); then echo 'Usage: ./remake.sh' >&2; exit 2; fi
if ! mkdir .run.lock 2>/dev/null; then
  echo 'A run is active or a stale .run.lock needs inspection; nothing removed.' >&2; exit 2
fi
trap 'rmdir .run.lock 2>/dev/null || true' EXIT
if [[ -e .archive.lock ]]; then
  echo 'An archive operation is active; nothing removed.' >&2; exit 2
fi
shopt -s nullglob
files=(app/*.x app/*.o app/*.mod app/*.smod app/*.a)
if (( ${#files[@]} )); then rm -f -- "${files[@]}"; fi
echo 'app/ compilation products removed. Run make to compile again.'
