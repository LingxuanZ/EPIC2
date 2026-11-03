# Run from the directory containing config.yml.
suppressPackageStartupMessages(library(Matrix))
config <- EPIC2::read_EPIC2_config("config.yml", check_files = TRUE)
set.seed(config$parameters$seed)

EPIC2::run_accessibility_stage(config, traits = "TC")

snp_accessibility <- readRDS("output/02_accessibility/TC.snp_by_cell.rds")
peak_accessibility <- readRDS("output/02_accessibility/TC.peak_by_cell.rds")
data.frame(
  branch = c("SNP", "Peak"),
  features = c(nrow(snp_accessibility), nrow(peak_accessibility)),
  cells = c(ncol(snp_accessibility), ncol(peak_accessibility))
)
