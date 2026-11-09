#!/usr/bin/env Rscript
# Usage: Rscript run_stages_2_to_5_EPIC2.R config.yml TC
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply the YAML path and one trait.")
trait <- args[2]
pipeline_started <- proc.time()[["elapsed"]]

log_event <- function(stage, status, elapsed = NULL, detail = NULL) {
  time <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  parts <- c(sprintf("[%s] [%s] [%s] %s", time, trait, stage, status),
             if (!is.null(elapsed)) sprintf("elapsed=%.1fs", elapsed),
             if (!is.null(detail)) gsub("[\r\n\t]+", " ", detail))
  cat(paste(parts, collapse = " "), "\n", sep = "")
  flush(stdout())
}

# Log each call without changing its arguments or suppressing package output.
run_step <- function(stage, expr) {
  started <- proc.time()[["elapsed"]]
  log_event(stage, "START")
  tryCatch({
    force(expr)
    NULL
  }, error = function(e) {
    log_event(stage, "FAILED", proc.time()[["elapsed"]] - started, conditionMessage(e))
    stop(e)
  }, interrupt = function(e) {
    log_event(stage, "INTERRUPTED", proc.time()[["elapsed"]] - started, conditionMessage(e))
    stop(e)
  })
  log_event(stage, "DONE", proc.time()[["elapsed"]] - started)
  invisible(NULL)
}

run_step("Setup", {
  suppressPackageStartupMessages({
    library(Matrix)
    library(SeuratObject)
    library(Signac)
  })
  config <- EPIC2::read_EPIC2_config(args[1], check_files = TRUE)
  if (!trait %in% names(config$gwas$traits)) stop("Unknown trait: ", trait)

  # The source files stay in the upload directory. Detect changes since submission.
  checksums <- file.path(config$output_dir, "source_checksums.rds")
  if (file.exists(checksums)) {
    expected <- readRDS(checksums)
    if (!identical(tools::md5sum(names(expected)), expected)) stop("R/YAML changed after submission.")
  }
  config$output_dir <- file.path(config$output_dir, trait)
  config$project$output_dir <- config$output_dir
  config$gwas$traits <- config$gwas$traits[trait]
  if (dir.exists(config$output_dir) && length(list.files(config$output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Trait output directory is not empty: ", config$output_dir)
  }

  # Start from the existing post-QC object; Stage 1 is not rerun.
  atac_path <- normalizePath(config$atac$input, mustWork = TRUE)
  atac <- readRDS(atac_path)
  if (!all(config$model$groupings %in% names(atac[[]]))) stop("Missing annotation columns in ATAC metadata.")
  fragments <- Signac::Fragments(atac[[config$atac$assay]])
  if (!length(fragments)) stop("The QC ATAC object must already have fragments attached.")
  fragment_paths <- vapply(fragments, Signac::GetFragmentData, character(1), slot = "path")
  if (any(file.access(c(fragment_paths, paste0(fragment_paths, ".tbi")), 4L) != 0L)) {
    stop("Attached fragment files or .tbi indexes are not readable on this compute node.")
  }
  rm(atac, fragments)
  invisible(gc())
  atac_dir <- file.path(config$output_dir, "01_atac")
  dir.create(atac_dir, recursive = TRUE)
  if (!file.symlink(atac_path, file.path(atac_dir, "atac.qc.rds"))) stop("Cannot link the QC ATAC input.")
  saveRDS(config, file.path(config$output_dir, "config.used.rds"))
  set.seed(config$parameters$seed)
})

# Within a trait, every call finishes before the next one starts.
run_step("Stage 2: GWAS", EPIC2::run_gwas_stage(config, traits = trait))
invisible(gc())
run_step("Stage 2: Accessibility", EPIC2::run_accessibility_stage(config, traits = trait))
invisible(gc())
run_step("Stage 3: LD", EPIC2::run_ld_stage(config, traits = trait))
invisible(gc())
run_step("Stage 4: WRS", EPIC2::run_wrs_stage(config, traits = trait))
invisible(gc())
run_step("Stage 5: Model", EPIC2::run_model_stage(config, traits = trait))

run_step("Finalize", {
  writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "sessionInfo.txt"))
  writeLines("COMPLETED", file.path(config$output_dir, "COMPLETED.txt"))
})
log_event("Pipeline", "DONE", proc.time()[["elapsed"]] - pipeline_started,
          paste0("output=", config$output_dir))
