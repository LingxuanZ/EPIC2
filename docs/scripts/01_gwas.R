# Run from the directory containing config.yml.
suppressPackageStartupMessages(library(Matrix))
config <- EPIC2::read_EPIC2_config("config.yml", check_files = TRUE)
set.seed(config$parameters$seed)

EPIC2::run_gwas_stage(config, gwas_name = "TC")

gwas_qc <- read.delim("output/02_gwas/TC.gwas_qc.tsv")
gwas_qc[, c("n_input", "n_reference_matched", "n_after_maf", "n_final")]
