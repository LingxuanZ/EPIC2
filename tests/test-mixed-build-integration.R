library(EPIC2)

# Mixed assemblies must fail before file reads or output-directory creation.
config <- EPIC2:::EPIC2_defaults()
config$genome$build <- "hg38"
config$atac$genome <- "hg19"
config$atac$input <- "not-read.rds"
config$gwas$traits <- list(TEST = list(summary = "not-read.tsv", genotype = "not-read.rds"))
config$output_dir <- tempfile("EPIC2-mixed-build-must-not-exist-")
message <- tryCatch({
  run_EPIC2(config)
  NA_character_
}, error = conditionMessage)
stopifnot(!is.na(message), grepl("atac.genome", message, fixed = TRUE), !dir.exists(config$output_dir))
message("Mixed-build rejection test passed.")
