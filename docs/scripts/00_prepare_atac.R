# Run from the directory containing config.yml.
config <- EPIC2::read_EPIC2_config("config.yml", check_files = TRUE)
set.seed(config$parameters$seed)

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
atac <- readRDS(config$atac$input)
fragment <- Signac::CreateFragmentObject(
  path = normalizePath("inputs/data/fragments.tsv.gz"),
  cells = setNames(colnames(atac), colnames(atac)),
  validate.fragments = FALSE, verbose = FALSE
)
Signac::Fragments(atac[["ATAC"]]) <- NULL
Signac::Fragments(atac[["ATAC"]]) <- fragment
checkpoint_dir <- file.path(config$output_dir, "01_atac")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(atac, file.path(checkpoint_dir, "atac.qc.rds"))
