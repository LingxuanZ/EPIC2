library(EPIC2)

run_common_build_integration <- function(genome_build = "hg19") {
  genome_build <- EPIC2:::assert_common_genome_build(genome_build)
  required <- c(
    "BiocParallel", "GenomeInfoDb", "GenomicRanges", "IRanges", "Matrix",
    "ggplot2", "presto", "Rsamtools", "Seurat", "Signac", "S4Vectors"
  )
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Common-build integration test is missing: ", paste(missing, collapse = ", "))
  }
  # Fail the test if the pipeline attempts to call the standalone converter.
  trace("lift_genomic_ranges", tracer = quote(stop("Unexpected pipeline liftOver call")),
        where = asNamespace("EPIC2"), print = FALSE)
  on.exit(untrace("lift_genomic_ranges", where = asNamespace("EPIC2")), add = TRUE)

  root_override <- Sys.getenv("EPIC2_INTEGRATION_OUTPUT")
  root <- if (nzchar(root_override)) root_override else tempfile("EPIC2-common-build-")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  if (!nzchar(root_override)) on.exit(unlink(root, recursive = TRUE), add = TRUE)
  cells <- paste0("cell", seq_len(8L))
  snp_position <- 1000000L + (0:5) * 10000L
  peaks <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = snp_position - 100L, end = snp_position + 100L)
  )
  GenomeInfoDb::genome(peaks) <- genome_build
  peak_ids <- paste0("chr1-", GenomicRanges::start(peaks), "-", GenomicRanges::end(peaks))
  names(peaks) <- peak_ids

  counts <- outer(seq_along(peaks), seq_along(cells), function(peak, cell) {
    1 + ((peak + cell + as.integer(cell > 4L) * peak) %% 3L)
  })
  rownames(counts) <- peak_ids
  colnames(counts) <- cells
  counts <- Matrix::Matrix(counts, sparse = TRUE)

  # FeatureMatrix counts fragment records, not the fifth-column read count.
  # Emit the intended number of distinct fragments for each cell and locus so
  # SNP coverage varies and all four association components are exercised.
  fragment_rows <- do.call(rbind, lapply(seq_along(peaks), function(peak) {
    do.call(rbind, lapply(seq_along(cells), function(cell) {
      offsets <- seq_len(as.integer(counts[peak, cell])) - 1L
      data.frame(
        chr = "chr1",
        start = snp_position[[peak]] - 50L + offsets,
        end = snp_position[[peak]] + 50L + offsets,
        cell = cells[[cell]],
        count = 1L,
        stringsAsFactors = FALSE
      )
    }))
  }))
  fragment_rows <- fragment_rows[order(fragment_rows$chr, fragment_rows$start), ]
  fragment_path <- file.path(root, "fragments.tsv")
  utils::write.table(
    fragment_rows, fragment_path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )
  fragment_gz <- Rsamtools::bgzip(fragment_path, overwrite = TRUE)
  Rsamtools::indexTabix(fragment_gz, format = "bed")

  chromatin <- Signac::CreateChromatinAssay(
    counts = counts,
    ranges = peaks,
    fragments = fragment_gz,
    validate.fragments = TRUE
  )
  metadata <- data.frame(
    tissue = rep(c("tissueA", "tissueB"), each = 4L),
    tissue_cell_type = rep(c("tissueA_type1", "tissueB_type1"), each = 4L),
    row.names = cells,
    stringsAsFactors = FALSE
  )
  atac <- Seurat::CreateSeuratObject(counts = chromatin, assay = "ATAC", meta.data = metadata)
  umap <- cbind(
    UMAP_1 = c(-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2),
    UMAP_2 = c(-1, -0.4, 0.5, 1, 1, 0.5, -0.4, -1)
  )
  rownames(umap) <- cells
  atac[["umap"]] <- Seurat::CreateDimReducObject(
    embeddings = umap, key = "UMAP_", assay = "ATAC"
  )
  atac$TSS.enrichment <- 5
  atac$FRIP <- 0.8
  atac$BlacklistRatio <- 0
  atac$nucleosome_signal <- 1
  atac_path <- file.path(root, "atac.rds")
  saveRDS(atac, atac_path)

  gwas <- data.frame(
    rsid = paste0("rs", seq_along(snp_position)),
    chr = 1L,
    pos = snp_position,
    A1 = rep(c("A", "C"), length.out = length(snp_position)),
    A2 = rep(c("G", "T"), length.out = length(snp_position)),
    MAF = 0.25,
    Beta = c(0.3, -0.2, 0.4, -0.1, 0.25, -0.35),
    SE = 0.1,
    P = c(0.01, 0.04, 0.001, 0.2, 0.02, 0.005),
    build = if (genome_build == "hg19") "GRCh37" else "GRCh38",
    stringsAsFactors = FALSE
  )
  genotype <- rbind(
    c(0, 0, 1, 1, 1, 2, 0, 1, 2, 0, 1, 2),
    c(0, 1, 0, 1, 2, 1, 2, 0, 1, 2, 1, 0),
    c(2, 1, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2),
    c(0, 0, 1, 2, 1, 0, 2, 1, 0, 1, 2, 1),
    c(1, 0, 2, 1, 0, 2, 1, 2, 0, 1, 0, 2),
    c(2, 2, 1, 0, 1, 0, 1, 0, 2, 1, 0, 1)
  )
  rownames(genotype) <- gwas$rsid
  colnames(genotype) <- paste0("reference", seq_len(ncol(genotype)))
  gwas_path <- file.path(root, "trait.tsv")
  genotype_path <- file.path(root, "genotype.tsv")
  utils::write.table(gwas, gwas_path, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(
    data.frame(rsid = rownames(genotype), genotype, check.names = FALSE),
    genotype_path, sep = "\t", quote = FALSE, row.names = FALSE
  )
  sldsc_path <- file.path(root, "sldsc_summary.tsv")
  utils::write.table(
    data.frame(
      trait = "SYNTHETIC",
      annotation_tested = c("tissueA_type1", "tissueB_type1"),
      coef_z = c(2.4, -1.3),
      Coefficient_P_value = c(0.016, 0.19),
      Enrichment = c(1.8, 0.7),
      Enrichment_p = c(0.025, 0.31),
      stringsAsFactors = FALSE
    ),
    sldsc_path, sep = "\t", quote = FALSE, row.names = FALSE
  )

  config <- EPIC2:::EPIC2_defaults()
  config$.config_dir <- root
  config$output_dir <- file.path(root, "results")
  config$genome$build <- genome_build
  config$atac$input <- atac_path
  config$atac$qc$tss_min <- 0
  config$atac$qc$frip_min <- 0
  config$atac$qc$count_quantile_min <- 0
  config$atac$qc$blacklist_max <- 1
  config$atac$qc$nucleosome_max <- 10
  config$atac$qc$min_cells_per_group <- 2L
  config$atac$qc$peak_width_min <- 50L
  config$atac$qc$peak_width_max <- 500L
  config$atac$qc$peak_remove_blacklist <- FALSE
  config$gwas$genome_column <- "build"
  config$gwas$maf_min <- 0.01
  config$gwas$ld_prune <- FALSE
  config$gwas$genotype_format <- "dosage"
  config$gwas$genotype_variant_id_column <- "rsid"
  config$gwas$traits <- list(
    SYNTHETIC = list(summary = gwas_path, genotype = genotype_path)
  )
  config$wrs$balanced$enabled <- TRUE
  config$wrs$balanced$cells_per_side <- 2L
  config$wrs$balanced$repeats <- 3L
  config$wrs$balanced$top_n <- 2L
  config$parameters$workers <- 1L
  config$model$groupings <- c("tissue_cell_type", "tissue")
  config$plots$sldsc_summary <- sldsc_path
  validate_EPIC2_config(config, check_files = TRUE)
  run_EPIC2(config)

  harmonized <- readRDS(file.path(config$output_dir, "02_gwas", "SYNTHETIC.gwas.rds"))
  peak_matrix <- readRDS(file.path(
    config$output_dir, "02_accessibility", "SYNTHETIC.peak_by_cell.rds"
  ))
  stopifnot(
    all(harmonized$EPIC2_input_genome == genome_build),
    all(harmonized$EPIC2_target_genome == genome_build),
    all(harmonized$pos == harmonized$EPIC2_input_pos),
    all(grepl("^chr1-", rownames(peak_matrix))),
    file.exists(file.path(
      config$output_dir, "05_models", "tissue_cell_type.cauchy_results.tsv"
    )),
    file.exists(file.path(config$output_dir, "figures", "figure_manifest.tsv")),
    file.exists(file.path(
      config$output_dir, "figures", "all_traits_GWAS_SNP_count_FeaturePlot.png"
    )),
    file.exists(file.path(
      config$output_dir, "figures", "all_traits.top3_per_group_topPeak_wilcoxon_DotPlot.pdf"
    )),
    file.exists(file.path(
      config$output_dir, "figures", "cauchy_vs_sldsc_aligned_heatmap_tissue_colored.png"
    )),
    file.exists(file.path(
      config$output_dir, "figures", "04_wrs", "balanced_wrs_pairwise_jaccard.png"
    )),
    file.exists(file.path(
      config$output_dir, "figures", "00_atac_qc", "peak_coverage_distribution.pdf"
    )),
    file.exists(file.path(
      config$output_dir, "figures", "02_accessibility",
      "all_traits_GWAS_SNP_count_VlnPlot.png"
    ))
  )
  figure_manifest <- utils::read.delim(
    file.path(config$output_dir, "figures", "figure_manifest.tsv"),
    stringsAsFactors = FALSE
  )
  stopifnot(
    !any(figure_manifest$status == "error"),
    !any(figure_manifest$status == "skipped")
  )
  for (grouping in config$model$groupings) {
    combined <- utils::read.delim(file.path(
      config$output_dir, "05_models", paste0(grouping, ".cauchy_results.tsv")
    ))
    component_columns <- c("p_snp_accessibility", "p_snp_wrs",
                           "p_peak_accessibility", "p_peak_wrs")
    component_p <- as.matrix(combined[, component_columns])
    stopifnot(
      nrow(combined) == 2L, all(is.finite(component_p)),
      all(component_p >= 0 & component_p <= 1),
      all(combined$n_tests_combined == 4L),
      all(is.finite(combined$p_cauchy)),
      all(combined$p_cauchy >= 0 & combined$p_cauchy <= 1)
    )
  }
  # Export switches must leave numerical outputs and plotting unchanged. Capture
  # ordinary/balanced WRS, complete LD models and final model results, excluding
  # receipts whose input fingerprints can change after a rerun.
  result_paths <- unlist(lapply(c("03_ld", "04_wrs", "05_models"), function(stage) {
    paths <- list.files(file.path(config$output_dir, stage), pattern = "\\.rds$", full.names = TRUE)
    paths[!grepl("(receipt|peak_response|peak_covariance)\\.rds$", paths)]
  }), use.names = FALSE)
  original_results <- lapply(result_paths, readRDS)
  wrs_rds <- list.files(file.path(config$output_dir, "04_wrs"),
                        pattern = "\\.rds$", full.names = TRUE)
  ordinary_wrs <- wrs_rds[!grepl("balanced|receipt", basename(wrs_rds))]
  optional_paths <- c(
    sub("\\.rds$", ".tsv", ordinary_wrs),
    file.path(config$output_dir, "03_ld",
              paste0("SYNTHETIC", c(".peak_response.rds", ".peak_covariance.rds")))
  )
  stopifnot(all(file.exists(optional_paths)))
  balanced_paths <- list.files(file.path(config$output_dir, "04_wrs"),
                               pattern = "balanced.*\\.tsv$", full.names = TRUE)
  balanced_before <- tools::md5sum(balanced_paths)
  # Cleanup must not remove exports for another trait or an unprocessed grouping.
  untouched_paths <- c(
    file.path(config$output_dir, "03_ld", "OTHER.peak_covariance.rds"),
    file.path(config$output_dir, "04_wrs", "SYNTHETIC.snp.unprocessed.tsv")
  )
  for (path in untouched_paths) writeLines("keep", path)

  config$ld$save_peak_checkpoints <- FALSE
  config$wrs$write_tsv <- FALSE
  # A configuration-only change must not invalidate completed statistical work.
  run_model_stage(config)
  # Failed reruns must not clean old optional exports or leave a valid receipt.
  for (stage in c("ld", "wrs")) {
    stage_function <- get(paste0("run_", stage, "_stage"), asNamespace("EPIC2"))
    override <- if (stage == "ld") "build_peak_ld_model_impl" else "run_wrs_matrix"
    environment(stage_function) <- list2env(
      setNames(list(function(...) stop("Simulated export rerun failure")), override),
      parent = asNamespace("EPIC2")
    )
    failure <- tryCatch({ stage_function(config); NA_character_ }, error = conditionMessage)
    stopifnot(identical(failure, "Simulated export rerun failure"),
              all(file.exists(optional_paths)),
              !file.exists(file.path(config$output_dir,
                if (stage == "ld") "03_ld" else "04_wrs",
                paste0("SYNTHETIC.", stage, "_receipt.rds"))))
  }
  run_EPIC2(config, stages = c("ld", "wrs", "model", "plot"))
  stopifnot(!any(file.exists(optional_paths)),
            all(file.exists(untouched_paths)),
            identical(tools::md5sum(balanced_paths), balanced_before))
  for (i in seq_along(result_paths)) {
    stopifnot(isTRUE(all.equal(readRDS(result_paths[[i]]), original_results[[i]])))
  }
  new_manifest <- utils::read.delim(file.path(config$output_dir, "figures", "figure_manifest.tsv"))
  stopifnot(identical(new_manifest$status, figure_manifest$status))
  recorded <- readRDS(file.path(config$output_dir, "config.used.rds"))
  stopifnot(identical(recorded$ld$save_peak_checkpoints, FALSE),
            identical(recorded$wrs$write_tsv, FALSE))
  unlink(untouched_paths)

  # Direct stage callers with old config objects still export, including after
  # the same output directory was previously run with exports disabled.
  legacy <- config
  legacy$ld$save_peak_checkpoints <- NULL
  legacy$wrs$write_tsv <- NULL
  run_ld_stage(legacy)
  run_wrs_stage(legacy)
  stopifnot(all(file.exists(optional_paths)))

  # A later sequential trait run must preserve the earlier trait's results.
  original <- readRDS(file.path(config$output_dir, "05_models", "tissue.cauchy_results.rds"))
  config$gwas$traits$SECOND <- config$gwas$traits$SYNTHETIC
  run_EPIC2(config, stages = c("gwas", "accessibility", "ld", "wrs", "model"), gwas_name = "SECOND")
  combined <- readRDS(file.path(config$output_dir, "05_models", "tissue.cauchy_results.rds"))
  stopifnot(setequal(unique(combined$trait), c("SYNTHETIC", "SECOND")),
            isTRUE(all.equal(combined$p_cauchy[combined$trait == "SYNTHETIC"], original$p_cauchy)))
  expect_failure <- function(pattern) {
    message <- tryCatch({ run_model_stage(config, traits = "SYNTHETIC"); NA_character_ }, error = conditionMessage)
    stopifnot(!is.na(message), grepl(pattern, message, fixed = TRUE))
  }
  runtime_path <- file.path(config$output_dir, "03_ld", "SYNTHETIC.ld_runtime.tsv")
  original_runtime <- readLines(runtime_path)
  writeLines(gsub("COMPLETED", "RUNNING", original_runtime, fixed = TRUE), runtime_path)
  expect_failure("LD is incomplete")
  writeLines(original_runtime, runtime_path)
  snp_path <- file.path(config$output_dir, "02_accessibility", "SYNTHETIC.snp_by_cell.rds")
  original_bytes <- readBin(snp_path, "raw", n = file.info(snp_path)$size)
  changed <- readRDS(snp_path); changed[1, 1] <- changed[1, 1] + 1
  saveRDS(changed, snp_path)
  expect_failure("changed")
  writeBin(original_bytes, snp_path)
  message("Common-build pipeline, plots, sequential traits, and incomplete/stale input guards passed: ", config$output_dir)
  invisible(TRUE)
}

if (identical(tolower(Sys.getenv("EPIC2_RUN_COMMON_BUILD_TEST")), "true")) {
  run_common_build_integration("hg19")
  run_common_build_integration("hg38")
} else {
  message("Skipping optional common-build integration test.")
}
