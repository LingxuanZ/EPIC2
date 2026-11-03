# Run from the directory containing config.yml.
suppressPackageStartupMessages(library(Matrix))
config <- EPIC2::read_EPIC2_config("config.yml", check_files = TRUE)
set.seed(config$parameters$seed)

EPIC2::run_model_stage(config, traits = "TC")

results <- read.delim("output/05_models/tissue_cell_type.cauchy_results.tsv")
head(results[order(results$p_cauchy),
  c("group", "p_cauchy", "p_cauchy_bh_within_trait")], 6)
