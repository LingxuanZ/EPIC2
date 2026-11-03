# Run from the directory containing config.yml.
suppressPackageStartupMessages(library(Matrix))
config <- EPIC2::read_EPIC2_config("config.yml", check_files = TRUE)
set.seed(config$parameters$seed)

EPIC2::run_wrs_stage(config, traits = "TC")

snp_wrs <- readRDS("output/04_wrs/TC.snp.tissue_cell_type.rds")
head(snp_wrs[, c("feature", "group", "pval_one_sided_enriched",
                 "z_one_sided_enriched_trunc")], 3)
