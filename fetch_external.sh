#!/usr/bin/env bash
# Fetch the unchanged upstream COUL90 dependency; it is not distributed here.
set -euo pipefail

usage() {
  echo 'Usage: ./fetch_external.sh [--source FILE]'
  echo 'Download COUL90 from its recorded upstream URL, or install an offline copy.'
  echo 'Every installed copy must match the fixed SHA-256 checksum.'
}
fail() { printf 'COUL90: %s\n' "$*" >&2; exit 1; }
source_file=''
case ${1:-} in
  '') (( $# == 0 )) || { usage >&2; exit 2; } ;;
  --help|-h) (( $# == 1 )) || { usage >&2; exit 2; }; usage; exit 0 ;;
  --source) (( $# == 2 )) || { usage >&2; exit 2; }; [[ -n "$2" ]] || fail 'Offline source path is empty.'; source_file=$2 ;;
  *) usage >&2; exit 2 ;;
esac
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
destination="$root/app/coulomb/coul90.f"
upstream="$root/app/coulomb/UPSTREAM.json"
expected='50ad795901114080ded8f63ebe3895b8460865f10fd3e59ff39d60e48456d13f'
command -v sha256sum >/dev/null 2>&1 || fail 'sha256sum is required.'
[[ -f "$upstream" ]] || fail "Missing source record: $upstream"
# These two scalar fields use the one-field-per-line JSON format kept in the repo.
source_url=$(awk -F '"' '$2 == "source_url" { print $4 }' "$upstream")
recorded_hash=$(awk -F '"' '$2 == "sha256" { print $4 }' "$upstream")
[[ "$recorded_hash" == "$expected" ]] || fail 'UPSTREAM.json does not match the pinned checksum.'
[[ "$source_url" == https://* && "$source_url" != *$'\n'* ]] || fail 'UPSTREAM.json must provide one HTTPS source_url.'
[[ ! -L "$destination" ]] || fail "Refusing to replace a symbolic link: $destination"
matching_hash() {
  local actual
  [[ -f "$1" ]] || return 1
  actual=$(sha256sum -- "$1") || return 1
  [[ "${actual%% *}" == "$expected" ]]
}
if matching_hash "$destination"; then
  echo 'COUL90 already installed; SHA-256 verified.'
  exit 0
fi
[[ ! -e "$destination" || -f "$destination" ]] || fail "Destination is not a regular file: $destination"
temporary=$(mktemp "$root/app/coulomb/.coul90.download.XXXXXX")
trap 'rm -f -- "$temporary"' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
if [[ -n "$source_file" ]]; then
  [[ -f "$source_file" ]] || fail "Offline source not found: $source_file"
  cp -- "$source_file" "$temporary" || fail 'Could not read the offline source; existing file preserved.'
else
  printf 'Downloading COUL90 from %s\n' "$source_url"
  if command -v curl >/dev/null 2>&1; then
    curl --disable --fail --location --proto '=https' --proto-redir '=https' \
      --connect-timeout 15 --max-time 90 --retry 2 --retry-delay 1 \
      --silent --show-error --output "$temporary" "$source_url" \
      || fail 'Download failed (see curl message); existing file preserved. Use --source FILE for an offline copy.'
  elif command -v wget >/dev/null 2>&1; then
    wget --no-config --no-hsts --https-only --timeout=30 --tries=2 --output-document="$temporary" "$source_url" \
      || fail 'Download failed (see wget message); existing file preserved. Use --source FILE for an offline copy.'
  else
    fail 'Install curl or wget, or supply an offline copy with --source FILE.'
  fi
fi
if ! matching_hash "$temporary"; then
  actual=$(sha256sum -- "$temporary")
  printf 'COUL90 checksum mismatch: expected %s; received %s\n' "$expected" "${actual%% *}" >&2
  fail 'Installation refused; existing file preserved.'
fi
chmod 644 -- "$temporary"
mv -f -- "$temporary" "$destination"
echo 'COUL90 installed; fixed SHA-256 verified; upstream bytes unchanged.'
