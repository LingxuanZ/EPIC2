library(EPIC2)

expect_error <- function(expression, pattern) {
  result <- tryCatch({ force(expression); NA_character_ }, error = conditionMessage)
  stopifnot(!is.na(result), grepl(pattern, result, fixed = TRUE))
}

# Every exported stage must reject an absent/conflicting common assembly before
# touching input files or creating its output directory.
for (stage in c("atac", "gwas", "accessibility", "ld", "wrs", "model", "plot")) {
  runner <- getExportedValue("EPIC2", paste0("run_", stage, "_stage"))
  config <- EPIC2:::EPIC2_defaults()
  config$atac$input <- "must-not-be-read.rds"
  config$gwas$traits <- list(TEST = list(summary = "unused.tsv", genotype = "unused.rds"))
  config$output_dir <- tempfile("EPIC2-stage-build-guard-")
  expect_error(runner(config), "genome.build")
  stopifnot(!dir.exists(config$output_dir))
  config$genome$build <- "hg38"
  config$atac$genome <- "hg19"
  expect_error(runner(config), "atac.genome")
  stopifnot(!dir.exists(config$output_dir))
}

# New matrix metadata prevents reuse under a conflicting assembly declaration.
matrix <- matrix(1, 2L, 2L)
attr(matrix, "genome_build") <- "hg19"
expect_error(EPIC2:::check_accessibility_matrix_build(matrix, "hg38"), "matrix_genome")
EPIC2:::check_accessibility_matrix_build(matrix, "GRCh37")

# A later-stage run must also reject a known mismatch in saved GWAS metadata,
# even when legacy configuration fields do not themselves conflict.
root <- tempfile("EPIC2-reused-build-")
dir.create(file.path(root, "02_gwas"), recursive = TRUE)
gwas <- data.frame(rsid = "rs1", chr = 1L, pos = 100L,
                   coord_id = "chr1-100-100", EPIC2_target_genome = "hg19")
saveRDS(gwas, file.path(root, "02_gwas", "TEST.gwas.rds"))
config <- EPIC2:::EPIC2_defaults()
config$genome$build <- "hg38"
config$gwas$traits <- list(TEST = list())
config$output_dir <- root
expect_error(run_ld_stage(config), "EPIC2_target_genome")
stopifnot(!file.exists(file.path(root, "03_ld", "TEST.ld_runtime.tsv")))
unlink(root, recursive = TRUE)
message("Common-build direct-stage and reused-input guards passed.")
