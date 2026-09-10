#!/usr/bin/env bash

# Run only stage 03 with the thread environment used by the Rmd LD job.
# Load the appropriate R module and set R_LIBS_USER before calling this script.
# Install the updated package separately; this script never installs packages.
#
# Usage:
#   bash scripts/run_ld_stage.sh "/path with spaces/config.yml" [HDL]
#   LD_THREADS=4 bash scripts/run_ld_stage.sh "/path/config.yml" HDL
#
# In Slurm, the default is SLURM_CPUS_PER_TASK; outside Slurm it is one thread.
# Set LD_THREADS explicitly to use a different number. It is capped at the
# Slurm allocation when SLURM_CPUS_PER_TASK is present.

set -euo pipefail

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  printf 'Usage: bash %s CONFIG.yml [TRAIT]\n' "$0"
  printf 'Default trait: HDL. Threads: LD_THREADS, then SLURM_CPUS_PER_TASK, then 1.\n'
  exit 0
fi

if (( $# < 1 || $# > 2 )); then
  printf 'Usage: bash %s CONFIG.yml [TRAIT]\n' "$0" >&2
  exit 2
fi

ld_config="$1"
ld_trait="${2:-HDL}"

if [[ ! -f "$ld_config" || ! -r "$ld_config" ]]; then
  printf 'Configuration file is missing or unreadable: %s\n' "$ld_config" >&2
  exit 2
fi

if [[ -z "$ld_trait" ]]; then
  printf 'The trait name must not be empty.\n' >&2
  exit 2
fi

ld_allocated_threads="${SLURM_CPUS_PER_TASK:-}"
ld_requested_threads="${LD_THREADS:-${ld_allocated_threads:-1}}"

if [[ ! "$ld_requested_threads" =~ ^[1-9][0-9]*$ ]]; then
  printf 'LD_THREADS must be a positive integer: %s\n' "$ld_requested_threads" >&2
  exit 2
fi

if [[ -n "$ld_allocated_threads" && ! "$ld_allocated_threads" =~ ^[1-9][0-9]*$ ]]; then
  printf 'SLURM_CPUS_PER_TASK must be a positive integer: %s\n' "$ld_allocated_threads" >&2
  exit 2
fi

ld_threads="$ld_requested_threads"

if [[ -n "$ld_allocated_threads" ]] && (( ld_threads > ld_allocated_threads )); then
  printf 'Limiting LD_THREADS=%s to SLURM_CPUS_PER_TASK=%s.\n' \
    "$ld_threads" "$ld_allocated_threads" >&2
  ld_threads="$ld_allocated_threads"
fi

if ! command -v Rscript >/dev/null 2>&1; then
  printf 'Rscript was not found. Load the same R module used for the Rmd first.\n' >&2
  exit 127
fi

# Set these BEFORE R starts: changing them inside an existing R session may
# not reconfigure an already initialized BLAS library. The LD loops remain
# serial; the linked BLAS/OpenMP implementation decides which operations
# actually use these threads. This is not four parallel R workers.
export OMP_NUM_THREADS="$ld_threads"
export MKL_NUM_THREADS="$ld_threads"
export OPENBLAS_NUM_THREADS="$ld_threads"

printf 'LD-only run started: %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
printf 'Configuration: %s\nTrait: %s\n' "$ld_config" "$ld_trait"
printf 'Slurm job: %s\nAllocated CPUs: %s\n' \
  "${SLURM_JOB_ID:-not set}" "${ld_allocated_threads:-not set}"
printf 'OMP_NUM_THREADS=%s MKL_NUM_THREADS=%s OPENBLAS_NUM_THREADS=%s\n' \
  "$OMP_NUM_THREADS" "$MKL_NUM_THREADS" "$OPENBLAS_NUM_THREADS"
printf 'R_LIBS_USER: %s\n' "${R_LIBS_USER:-not set; using R defaults}"

exec Rscript --vanilla - "$ld_config" "$ld_trait" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)

if (!requireNamespace("EPIC2", quietly = TRUE)) {
  stop(
    "EPIC2 cannot be loaded from the selected R libraries. ",
    "Install the updated local EPIC2 source into R_LIBS_USER before submission.\n",
    "Libraries searched: ", paste(.libPaths(), collapse = "; "),
    call. = FALSE
  )
}

message("R: ", R.version.string)
message("EPIC2 version: ", as.character(utils::packageVersion("EPIC2")))
message("EPIC2 library: ", find.package("EPIC2"))
message("BLAS: ", unname(extSoftVersion()["BLAS"]))

# Reuse the existing configuration and completed GWAS/accessibility outputs.
# Do not rerun upstream stages or change statistical LD parameters.
config <- EPIC2::read_EPIC2_config(args[[1L]])
EPIC2::run_EPIC2(config, stages = "ld", gwas_name = args[[2L]])

message("LD-only run finished: ", format(Sys.time(), tz = "UTC", usetz = TRUE))
RSCRIPT
