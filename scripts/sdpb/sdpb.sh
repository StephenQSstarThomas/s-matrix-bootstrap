#!/usr/bin/env bash
# Run an SDPB executable (sdpb / pmp2sdp / spectrum / approx_objective) in the
# local container against a working directory.
#
#   scripts/sdpb/sdpb.sh WORKDIR pmp2sdp --precision=768 --input=pmp.json --output=sdp
#   scripts/sdpb/sdpb.sh WORKDIR sdpb    --precision=768 -s sdp -o out
#   scripts/sdpb/sdpb.sh WORKDIR -n 16 sdpb --precision=768 -s sdp -o out
#
# WORKDIR is bind-mounted at /w and is the container's cwd, so every path you
# pass to the executable is relative to WORKDIR.
#
# Do NOT add --network none: SDPB initialises MPI at startup and hangs forever
# with no network namespace.  The daemon's iptables rules already confine it.
set -euo pipefail

# SDPB parallelizes over MPI ranks; its BLAS calls must be single-threaded.
# Keep this local to the launcher so Python operator assembly may use threads.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

image="${SDPB_IMAGE:-bootstrapcollaboration/sdpb:3.1.0}"
[[ $# -ge 2 ]] || { echo "usage: $0 WORKDIR [-n NPROC] EXECUTABLE [args...]" >&2; exit 2; }
workdir="$(realpath -e -- "$1")"; shift
[[ "$workdir" != *,* ]] || { echo "workdir cannot contain a comma" >&2; exit 2; }

nproc_args=()
if [[ "${1:-}" == -n ]]; then
  nproc_args=(mpirun -n "$2" --allow-run-as-root
    -x OPENBLAS_NUM_THREADS -x OMP_NUM_THREADS -x MKL_NUM_THREADS)
  shift 2
fi

# A real cluster uses an installed SDPB/MPI module on a shared filesystem.
# Docker's local -n counts ranks on ONE host; never invent remote hosts.
if [[ "${SDPB_NATIVE:-0}" == 1 ]]; then
  if [[ -n "${SDPB_HOSTFILE:-}" && ${#nproc_args[@]} -gt 0 ]]; then
    [[ -f "$SDPB_HOSTFILE" ]] || { echo "missing SDPB_HOSTFILE" >&2; exit 2; }
    nproc_args+=(--hostfile "$SDPB_HOSTFILE")
  fi
  cd "$workdir"
  exec "${nproc_args[@]}" "$@"
fi
[[ -z "${SDPB_HOSTFILE:-}" ]] || { echo "remote hosts require SDPB_NATIVE=1" >&2; exit 2; }
cid_args=()
if [[ -n "${SDPB_CIDFILE:-}" ]]; then cid_args=(--cidfile "$SDPB_CIDFILE"); fi

# SDPB allocates MPI shared-memory windows; Docker's default 64 MB /dev/shm
# makes MPI_Win_allocate_shared fail with MPI_ERR_INTERN once the problem is
# more than a toy.  Size it from SDPB_SHM (default 16g).
exec docker run --rm --user "$(id -u):$(id -g)" \
  "${cid_args[@]}" \
  --env OPENBLAS_NUM_THREADS=1 --env OMP_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --shm-size "${SDPB_SHM:-16g}" \
  --mount "type=bind,src=$workdir,dst=/w" --workdir /w \
  "$image" "${nproc_args[@]}" "$@"
