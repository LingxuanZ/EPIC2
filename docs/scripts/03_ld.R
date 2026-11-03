# Run from the directory containing config.yml.
suppressPackageStartupMessages(library(Matrix))
config <- EPIC2::read_EPIC2_config("config.yml", check_files = TRUE)
set.seed(config$parameters$seed)

EPIC2::run_ld_stage(config, traits = "TC")

ld_status <- read.delim("output/03_ld/TC.ld_runtime.tsv")
ld_status[ld_status$setting == "status", , drop = FALSE]
snp_blocks <- readRDS("output/03_ld/TC.snp_blocks.rds")
data.frame(
  SNP_blocks = length(snp_blocks),
  largest_block = max(lengths(snp_blocks))
)
