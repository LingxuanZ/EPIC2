#!/bin/bash
# Upload this file, run_stages_2_to_5_EPIC2.R and config.yml together.
# Run on the login node: bash submit_stages_2_to_5_EPIC2.sh
set -euo pipefail

EPIC2_SLURM_ACCOUNT="YOUR_SLURM_ACCOUNT"
if [[ -z "${EPIC2_SLURM_ACCOUNT}" || "${EPIC2_SLURM_ACCOUNT}" == "YOUR_SLURM_ACCOUNT" ]]; then
  echo "Set EPIC2_SLURM_ACCOUNT at the top of this script before running." >&2
  exit 1
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/config.yml"
WORKER="${SCRIPT_DIR}/run_stages_2_to_5_EPIC2.R"
EPIC2_REF="e0a1db5dabd261c582b90ad6b2ecb221d86d9492"
PACKAGE_LIBRARY="${HOME}/.local/share/EPIC2/${EPIC2_REF}/R_library"
BASE_LIBRARY="/path/to/R/library"
TRAITS=(HDL LDL TC TG SCZ SCZ2 BIP SCZBIP)
[[ -r "${CONFIG_FILE}" && -r "${WORKER}" ]] || { echo "Missing YAML or R script." >&2; exit 1; }

module purge
module load GCC/14.2.0 R/4.5.1
export R_LIBS_USER="${BASE_LIBRARY}"
export R_LIBS="${BASE_LIBRARY}"
RUN_DIR=$(Rscript --vanilla - "${CONFIG_FILE}" "${TRAITS[@]}" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
config <- yaml::read_yaml(args[1])
stopifnot(identical(names(config$gwas$traits), args[-1]))
path <- config$project$output_dir
if (length(path) != 1L || !startsWith(path, "/")) stop("Use an absolute project.output_dir in YAML.")
inputs <- c(config$atac$input, unlist(lapply(config$gwas$traits, function(x) c(x$summary, x$genotype))))
missing <- inputs[!file.exists(inputs) | file.access(inputs, 4L) != 0L]
if (length(missing)) stop("Missing or unreadable input(s): ", paste(missing, collapse = ", "))
cat(path)
RS
)

# Existing empty directories are allowed. Preserve prior submissions and outputs.
if [[ -e "${RUN_DIR}/array_job_id.txt" ]]; then
  echo "An array was already submitted: $(cat "${RUN_DIR}/array_job_id.txt"). Check it before reusing ${RUN_DIR}." >&2
  exit 1
fi
for trait in "${TRAITS[@]}"; do
  if [[ -e "${RUN_DIR}/${trait}" || -L "${RUN_DIR}/${trait}" ]]; then
    echo "Trait output already exists: ${RUN_DIR}/${trait}. Choose a fresh output root in YAML." >&2
    exit 1
  fi
done
mkdir -p "${RUN_DIR}/log" "${PACKAGE_LIBRARY}"
export R_LIBS_USER="${PACKAGE_LIBRARY}"
export R_LIBS="${PACKAGE_LIBRARY}:${BASE_LIBRARY}"

# Cache the pinned package under HOME; do not copy the R/YAML/shell files to scratch.
Rscript --vanilla - "${CONFIG_FILE}" "${WORKER}" "${PACKAGE_LIBRARY}" "${EPIC2_REF}" \
  <<'RS' 2>&1 | tee "${RUN_DIR}/log/installation.log"
args <- commandArgs(trailingOnly = TRUE)
installed <- suppressWarnings(utils::packageDescription("EPIC2", lib.loc = args[3]))
if (!is.list(installed) || !identical(installed$RemoteSha, args[4])) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  remotes::install_github(paste0("LingxuanZ/EPIC2@", args[4]), lib = args[3],
    dependencies = FALSE, upgrade = "never", build_vignettes = FALSE)
}
stopifnot(identical(utils::packageDescription("EPIC2")$RemoteSha, args[4]))
packages <- c("Matrix", "Seurat", "SeuratObject", "Signac", "GenomicRanges", "GenomeInfoDb",
              "IRanges", "S4Vectors", "BiocParallel", "Rsamtools", "igraph", "presto")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing existing dependencies: ", paste(missing, collapse = ", "))
config <- EPIC2::read_EPIC2_config(args[1], check_files = TRUE)
saveRDS(tools::md5sum(args[1:2]), file.path(config$output_dir, "source_checksums.rds"))
writeLines(args[4], file.path(config$output_dir, "EPIC2_COMMIT.txt"))
cat("EPIC2 commit:", args[4], "\nLibrary:", find.package("EPIC2"), "\n")
RS

export FETAL_PIPELINE_DIR="${SCRIPT_DIR}"
export FETAL_PIPELINE_LIBRARY="${PACKAGE_LIBRARY}"
export FETAL_PIPELINE_BASE_LIBRARY="${BASE_LIBRARY}"
cd -- "${RUN_DIR}"
SUBMISSION=$(/usr/bin/sbatch --account="${EPIC2_SLURM_ACCOUNT}" --parsable --export=ALL --chdir=. \
  --output=log/%x.%A_%a.out --error=log/%x.%A_%a.err <<'SBATCH'
#!/bin/bash
#SBATCH --job-name=epic2_stages_2_to_5
#SBATCH --array=1-8
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=160G
set -euo pipefail
TRAITS=(HDL LDL TC TG SCZ SCZ2 BIP SCZBIP)
[[ "${SLURM_ARRAY_TASK_ID:-}" =~ ^[1-8]$ ]] || { echo "Invalid array task ID." >&2; exit 1; }
TRAIT="${TRAITS[$((SLURM_ARRAY_TASK_ID - 1))]}"
EPIC2_TASK_STARTED="$(date +%s)"
log_task_exit() {
  local exit_code=$?
  local status="FAILED"
  [[ "${exit_code}" -eq 0 ]] && status="DONE"
  printf '[%s] [%s] [Task] %s elapsed=%ss exit=%s\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${TRAIT}" "${status}" \
    "$(( $(date +%s) - EPIC2_TASK_STARTED ))" "${exit_code}"
  exit "${exit_code}"
}
trap log_task_exit EXIT
printf '[%s] [%s] [Task] START job=%s array_task=%s\n' \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${TRAIT}" "${SLURM_JOB_ID:-unknown}" "${SLURM_ARRAY_TASK_ID}"

module purge
module load GCC/14.2.0 R/4.5.1
export R_LIBS_USER="${FETAL_PIPELINE_LIBRARY}"
export R_LIBS="${FETAL_PIPELINE_LIBRARY}:${FETAL_PIPELINE_BASE_LIBRARY}"
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
Rscript --vanilla "${FETAL_PIPELINE_DIR}/run_stages_2_to_5_EPIC2.R" \
  "${FETAL_PIPELINE_DIR}/config.yml" "${TRAIT}"
SBATCH
)
ARRAY_JOB="${SUBMISSION%%;*}"
[[ "${ARRAY_JOB}" =~ ^[0-9]+$ ]] || { echo "Unexpected Slurm response: ${SUBMISSION}" >&2; exit 1; }
printf '%s\n' "${ARRAY_JOB}" > "${RUN_DIR}/array_job_id.txt"
echo "Submitted stages 2-5 array: ${ARRAY_JOB} (eight traits)"
echo "Results: ${RUN_DIR}/<TRAIT>/"
echo "Model results: ${RUN_DIR}/<TRAIT>/05_models/"
echo "Monitor: squeue -j ${ARRAY_JOB}"
