library(EPIC2)

run_mixed_build_integration <- function() {
  required <- c(
    "BiocParallel", "GenomeInfoDb", "GenomicRanges", "IRanges", "Matrix",
    "ggplot2", "presto", "Rsamtools", "Seurat", "Signac", "S4Vectors", "rtracklayer"
  )
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Mixed-build integration test is missing: ", paste(missing, collapse = ", "))
  }
  chain <- Sys.getenv("EPIC2_HG19_TO_HG38_CHAIN")
  if (!nzchar(chain) || !file.exists(chain)) {
    stop("Set EPIC2_HG19_TO_HG38_CHAIN to a readable directional UCSC chain file.")
  }

  root_override <- Sys.getenv("EPIC2_INTEGRATION_OUTPUT")
  root <- if (nzchar(root_override)) root_override else tempfile("EPIC2-mixed-build-")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  if (!nzchar(root_override)) on.exit(unlink(root, recursive = TRUE), add = TRUE)
  cells <- paste0("cell", seq_len(8L))
  snp_position <- 1000000L + (0:5) * 10000L
  peaks <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = snp_position - 100L, end = snp_position + 100L)
  )
  GenomeInfoDb::genome(peaks) <- "hg19"
  peak_ids <- paste0("chr1-", GenomicRanges::start(peaks), "-", GenomicRanges::end(peaks))
  names(peaks) <- peak_ids

  counts <- outer(seq_along(peaks), seq_along(cells), function(peak, cell) {
    1 + ((peak + cell + as.integer(cell > 4L) * peak) %% 3L)
  })
  rownames(counts) <- peak_ids
  colnames(counts) <- cells
  counts <- Matrix::Matrix(counts, sparse = TRUE)

  fragment_rows <- do.call(rbind, lapply(seq_along(peaks), function(peak) {
    data.frame(
      chr = "chr1",
      start = snp_position[[peak]] - 50L,
      end = snp_position[[peak]] + 50L,
      cell = cells,
      count = as.integer(counts[peak, ]),
      stringsAsFactors = FALSE
    )
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
    build = "GRCh37",
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
  config$genome$target <- "hg38"
  config$genome$chains$hg19_to_hg38 <- chain
  config$atac$input <- atac_path
  config$atac$genome <- "hg19"
  config$atac$qc$tss_min <- 0
  config$atac$qc$frip_min <- 0
  config$atac$qc$count_quantile_min <- 0
  config$atac$qc$blacklist_max <- 1
  config$atac$qc$nucleosome_max <- 10
  config$atac$qc$min_cells_per_group <- 2L
  config$atac$qc$peak_width_min <- 50L
  config$atac$qc$peak_width_max <- 500L
  config$atac$qc$peak_remove_blacklist <- FALSE
  config$gwas$genome <- "hg19"
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
    all(harmonized$EPIC2_input_genome == "hg19"),
    all(harmonized$EPIC2_target_genome == "hg38"),
    all(harmonized$pos != harmonized$EPIC2_input_pos),
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
  message("Mixed-build pipeline and publication-figure integration test completed: ", config$output_dir)
  invisible(TRUE)
}

if (identical(tolower(Sys.getenv("EPIC2_RUN_MIXED_BUILD_TEST")), "true")) {
  run_mixed_build_integration()
} else {
  message("Skipping optional mixed-build integration test.")
}
