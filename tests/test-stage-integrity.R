library(EPIC2)

local({
  root <- tempfile("epic2-integrity-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  expect_error <- function(expr, pattern) {
    message <- tryCatch({ force(expr); NA_character_ }, error = conditionMessage)
    stopifnot(!is.na(message), grepl(pattern, message, fixed = TRUE))
  }
  input <- file.path(root, "input"); output <- file.path(root, "output")
  writeLines("input1", input); writeLines("output1", output)
  path <- file.path(root, "receipt.rds")
  settings <- list(ridge = 0.01)
  EPIC2:::write_stage_receipt(path, EPIC2:::file_fingerprints(input), output, settings)
  EPIC2:::validate_stage_receipt(path, settings)
  expect_error(EPIC2:::validate_stage_receipt(path, list(ridge = 0.02)), "settings changed")
  writeLines("input2", input)
  expect_error(EPIC2:::validate_stage_receipt(path, settings), "inputs changed")
  writeLines("input1", input); writeLines("output2", output)
  expect_error(EPIC2:::validate_stage_receipt(path, settings), "outputs changed")

  config <- EPIC2:::EPIC2_defaults()
  config$genome$build <- "hg19"
  config$output_dir <- file.path(root, "results")
  config$.config_dir <- root
  config$atac$input <- "atac.rds"
  config$gwas$traits <- list(TEST = list(summary = "gwas.tsv", genotype = "genotype.rds"))
  lock <- EPIC2:::acquire_stage_lock(config, "ld")
  expect_error(EPIC2:::acquire_stage_lock(config, "model"), "Another stage")
  unlink(lock, recursive = TRUE)
  config$ld$ridge <- 0.2
  EPIC2:::write_run_manifest(config, "ld", "TEST")
  restored <- read_EPIC2_config(file.path(config$output_dir, "config.used.yml"))
  stopifnot(restored$ld$ridge == 0.2,
            restored$atac$input == file.path(root, "atac.rds"),
            restored$gwas$traits$TEST$summary == file.path(root, "gwas.tsv"),
            normalizePath(restored$output_dir) == normalizePath(config$output_dir),
            length(list.files(file.path(config$output_dir, "run_history"))) == 2L)
  expect_error(validate_genome_coordinates(NA, 100, genome = "hg19"), "invalid")
  expect_error(validate_genome_coordinates("chr1", 100.5, genome = "hg19"), "invalid")
})
message("Completion checksums, stage exclusion, effective configuration, and coordinate validation passed.")
