# Run from the directory containing config.yml.
suppressPackageStartupMessages(library(Matrix))
config <- EPIC2::read_EPIC2_config("config.yml", check_files = TRUE)
set.seed(config$parameters$seed)

EPIC2::run_plot_stage(config, traits = "TC")

figure_manifest <- read.delim("output/figures/figure_manifest.tsv")
table(figure_manifest$status)
