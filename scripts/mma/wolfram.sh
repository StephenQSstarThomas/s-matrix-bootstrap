#!/usr/bin/env bash
# Run a Wolfram Language script in the licensed container.
#
# Role: independent checks of the published formulas using our Mathematica
# scripts, not a recovered author notebook. The SDPB CLI can call this audit
# before assembly; Mathematica does not solve the bootstrap optimization.
#
# Usage: scripts/mma/wolfram.sh SCRIPT.wls [ro:HOST:DST ...] [rw:HOST:DST ...] [-- ARGS...]
#        Everything after a literal -- is passed to the script as $ScriptCommandLine.
#   Env: MMA_LICENSE_DIR (default /home/shiqiu/Licensing), MMA_IMAGE.
set -euo pipefail

license_dir="${MMA_LICENSE_DIR:-/home/shiqiu/Licensing}"
image="${MMA_IMAGE:-wolframresearch/wolframengine:15.0.0}"
[[ -f "$license_dir/mathpass" ]] || { echo "no mathpass in $license_dir" >&2; exit 2; }
[[ $# -ge 1 ]] || { echo "usage: $0 SCRIPT.wls [ro:HOST:DST ...] [rw:HOST:DST ...]" >&2; exit 2; }

script="$(realpath -e -- "$1")"; shift
mounts=(
  --mount "type=bind,src=$license_dir,dst=/licenses,readonly"
  --mount "type=bind,src=$(dirname -- "$script"),dst=/work,readonly"
)
script_args=()
seen_sep=0
for spec in "$@"; do
  if [[ $seen_sep == 1 ]]; then script_args+=("$spec"); continue; fi
  if [[ "$spec" == -- ]]; then seen_sep=1; continue; fi
  mode="${spec%%:*}"; rest="${spec#*:}"; host="${rest%%:*}"; dst="${rest#*:}"
  host="$(realpath -e -- "$host")"
  # Docker --mount separates options with commas; a comma in a path is ambiguous.
  [[ "$host$dst" != *,* ]] || { echo "mount paths cannot contain commas" >&2; exit 2; }
  case "$mode" in
    ro) mounts+=(--mount "type=bind,src=$host,dst=$dst,readonly") ;;
    rw) mounts+=(--mount "type=bind,src=$host,dst=$dst") ;;
    *)  echo "mount spec must start with ro: or rw: -- got $spec" >&2; exit 2 ;;
  esac
done

# --network none keeps the container offline AND makes $MachineID stable, which
# is what lets the container-bound mathpass activate on every run.
# --user keeps files written to rw: mounts owned by the caller, not root.
cid_args=()
if [[ -n "${MMA_CIDFILE:-}" ]]; then cid_args=(--cidfile "$MMA_CIDFILE"); fi
exec docker run --rm -i --network none --user "$(id -u):$(id -g)" \
  "${cid_args[@]}" \
  "${mounts[@]}" \
  -e WOLFRAMINIT="-pwfile /licenses/mathpass" -e HOME=/tmp \
  --workdir /work \
  "$image" \
  /usr/bin/wolframscript -file "/work/$(basename -- "$script")" "${script_args[@]}"
