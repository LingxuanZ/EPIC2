make_bpparam <- function(workers) {
  require_namespaces("BiocParallel", "parallel accessibility calculation")
  workers <- max(1L, as.integer(workers))
  if (workers == 1L) return(BiocParallel::SerialParam())
  if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(workers = workers, type = "SOCK")
  } else {
    BiocParallel::MulticoreParam(workers = workers)
  }
}

gwas_ranges <- function(gwas, columns = EPIC2_defaults()$gwas$columns) {
  require_namespaces(c("GenomicRanges", "IRanges"), "GWAS genomic ranges")
  coord <- gwas$coord_id %||% coordinate_id(gwas[[columns$chr]], gwas[[columns$pos]])
  ranges <- GenomicRanges::GRanges(
    seqnames = normalize_chr(gwas[[columns$chr]]),
    ranges = IRanges::IRanges(
      start = as.integer(gwas[[columns$pos]]),
      end = as.integer(gwas[[columns$pos]])
    ),
    strand = "*"
  )
  names(ranges) <- coord
  ranges
}

prepare_gwas_ranges_for_atac <- function(gwas, columns, gwas_genome, atac_genome,
                                         chain_file = NULL,
                                         unmapped_action = "drop",
                                         multimapping_action = "drop",
                                         duplicate_action = "drop") {
  input_builds <- unique(as.character(gwas$EPIC2_input_genome %||% character()))
  input_builds <- input_builds[!is.na(input_builds) & nzchar(input_builds)]
  use_original <- all(c("EPIC2_input_chr", "EPIC2_input_pos") %in% colnames(gwas)) &&
    length(input_builds) == 1L && normalize_genome_build(input_builds) == normalize_genome_build(atac_genome)
  if (use_original) {
    ranges <- GenomicRanges::GRanges(
      seqnames = canonical_chr(gwas$EPIC2_input_chr),
      ranges = IRanges::IRanges(
        start = as.integer(gwas$EPIC2_input_pos),
        end = as.integer(gwas$EPIC2_input_pos)
      ),
      strand = "*"
    )
    names(ranges) <- gwas$coord_id
    return(lift_genomic_ranges(
      ranges, atac_genome, atac_genome,
      unmapped_action = unmapped_action,
      multimapping_action = multimapping_action,
      duplicate_action = duplicate_action,
      label = "original GWAS variants reused in the matching ATAC build"
    ))
  }
  ranges <- gwas_ranges(gwas, columns)
  lift_genomic_ranges(
    ranges,
    source_genome = gwas_genome,
    target_genome = atac_genome,
    chain_file = chain_file,
    unmapped_action = unmapped_action,
    multimapping_action = multimapping_action,
    duplicate_action = duplicate_action,
    label = "target-build GWAS variants mapped to the ATAC build"
  )
}

#' Build a trait-specific SNP-by-cell accessibility matrix
#'
#' @param atac Filtered Seurat/Signac object with fragments attached.
#' @param gwas Harmonized GWAS table containing `coord_id`.
#' @param columns GWAS column mapping.
#' @param assay ATAC assay name.
#' @param workers Number of BiocParallel workers.
#' @param chunks Number of genomic feature chunks.
#' @param gwas_genome Build used by the harmonized GWAS coordinates.
#' @param atac_genome Build used by ATAC peaks and fragments.
#' @param chain_file Directional chain from `gwas_genome` to `atac_genome`.
#' @param unmapped_action Drop or reject unmapped ranges.
#' @param multimapping_action Drop or reject multiply mapped ranges.
#' @param duplicate_action Drop or reject duplicate target ranges.
#' @return Sparse SNP-by-cell matrix containing nonzero SNPs in ATAC peaks.
#' @export
build_snp_accessibility <- function(atac, gwas, columns = EPIC2_defaults()$gwas$columns,
                                    assay = "ATAC", workers = 1L, chunks = 20L,
                                    gwas_genome = "hg19", atac_genome = gwas_genome,
                                    chain_file = NULL, unmapped_action = "drop",
                                    multimapping_action = "drop", duplicate_action = "drop") {
  require_namespaces(
    c("Signac", "GenomicRanges", "GenomeInfoDb", "BiocParallel", "Matrix"),
    "SNP accessibility"
  )
  counts <- get_atac_counts(atac, assay)
  peak_ranges <- Signac::StringToGRanges(rownames(counts), sep = c("-", "-"))
  peak_ranges <- normalize_granges_seqnames(peak_ranges)
  validate_genome_coordinates(
    as.character(GenomicRanges::seqnames(peak_ranges)), GenomicRanges::start(peak_ranges),
    GenomicRanges::end(peak_ranges), atac_genome, "ATAC peaks"
  )
  prepared <- prepare_gwas_ranges_for_atac(
    gwas, columns, gwas_genome, atac_genome, chain_file,
    unmapped_action, multimapping_action, duplicate_action
  )
  snps <- prepared$ranges
  snps <- IRanges::subsetByOverlaps(snps, peak_ranges, ignore.strand = TRUE)
  if (length(snps) == 0L) stop("No harmonized GWAS variants overlap retained ATAC peaks.", call. = FALSE)

  n_chunks <- min(max(1L, as.integer(chunks)), length(snps))
  split_id <- cut(seq_along(snps), breaks = n_chunks, labels = FALSE)
  feature_chunks <- split(snps, split_id)
  param <- make_bpparam(workers)
  matrices <- BiocParallel::bplapply(feature_chunks, function(features) {
    matrix <- Signac::FeatureMatrix(
      fragments = Signac::Fragments(atac[[assay]]),
      features = features,
      cells = colnames(atac)
    )
    rownames(matrix) <- names(features)
    matrix
  }, BPPARAM = param)
  matrix <- do.call(rbind, matrices)
  matrix <- matrix[Matrix::rowSums(matrix) > 0, , drop = FALSE]
  matrix <- matrix[gwas$coord_id[gwas$coord_id %in% rownames(matrix)], , drop = FALSE]
  if (nrow(matrix) == 0L) {
    stop("No SNP accessibility rows contain ATAC fragments after genome harmonization.", call. = FALSE)
  }
  attr(matrix, "genome_qc") <- prepared$qc
  matrix
}

#' Build a trait-specific peak-by-cell accessibility matrix
#'
#' @param atac Filtered Seurat/Signac object.
#' @param gwas Harmonized GWAS table.
#' @param columns GWAS column mapping.
#' @param assay ATAC assay name.
#' @param gwas_genome Build used by the harmonized GWAS coordinates.
#' @param atac_genome Build used by ATAC peaks and fragments.
#' @param gwas_to_atac_chain Directional chain used to query source-build fragments.
#' @param atac_to_gwas_chain Directional chain used to name retained peaks in the target build.
#' @param unmapped_action Drop or reject unmapped ranges.
#' @param multimapping_action Drop or reject multiply mapped ranges.
#' @param duplicate_action Drop or reject duplicate target ranges.
#' @return A list with a sparse peak matrix and SNP-to-peak overlap table.
#' @export
build_peak_accessibility <- function(atac, gwas, columns = EPIC2_defaults()$gwas$columns,
                                     assay = "ATAC", gwas_genome = "hg19",
                                     atac_genome = gwas_genome,
                                     gwas_to_atac_chain = NULL, atac_to_gwas_chain = NULL,
                                     unmapped_action = "drop",
                                     multimapping_action = "drop",
                                     duplicate_action = "drop") {
  require_namespaces(c("Signac", "GenomicRanges", "GenomeInfoDb"), "peak accessibility")
  counts <- get_atac_counts(atac, assay)
  peaks <- Signac::StringToGRanges(rownames(counts), sep = c("-", "-"))
  peaks <- normalize_granges_seqnames(peaks)
  names(peaks) <- rownames(counts)
  validate_genome_coordinates(
    as.character(GenomicRanges::seqnames(peaks)), GenomicRanges::start(peaks),
    GenomicRanges::end(peaks), atac_genome, "ATAC peaks"
  )
  prepared <- prepare_gwas_ranges_for_atac(
    gwas, columns, gwas_genome, atac_genome, gwas_to_atac_chain,
    unmapped_action, multimapping_action, duplicate_action
  )
  snps <- prepared$ranges
  hits <- GenomicRanges::findOverlaps(snps, peaks, ignore.strand = TRUE)
  if (length(hits) == 0L) stop("No harmonized GWAS variants overlap retained ATAC peaks.", call. = FALSE)
  source_overlap <- unique(data.frame(
    snp = names(snps)[S4Vectors::queryHits(hits)],
    peak = names(peaks)[S4Vectors::subjectHits(hits)],
    stringsAsFactors = FALSE
  ))
  source_peak_order <- rownames(counts)[rownames(counts) %in% unique(source_overlap$peak)]
  peak_lift <- lift_genomic_ranges(
    peaks[source_peak_order],
    source_genome = atac_genome,
    target_genome = gwas_genome,
    chain_file = atac_to_gwas_chain,
    unmapped_action = unmapped_action,
    multimapping_action = multimapping_action,
    duplicate_action = duplicate_action,
    label = "GWAS-overlapping ATAC peaks mapped to the target build"
  )
  retained_source_peaks <- names(peak_lift$ranges)
  target_peak_ids <- paste0(
    canonical_chr(as.character(GenomicRanges::seqnames(peak_lift$ranges))), "-",
    GenomicRanges::start(peak_lift$ranges), "-", GenomicRanges::end(peak_lift$ranges)
  )
  source_to_target <- setNames(target_peak_ids, retained_source_peaks)
  source_overlap <- source_overlap[source_overlap$peak %in% retained_source_peaks, , drop = FALSE]
  overlap <- source_overlap
  overlap$peak <- unname(source_to_target[overlap$peak])
  overlap <- unique(overlap)
  source_peak_order <- source_peak_order[source_peak_order %in% retained_source_peaks]
  matrix <- counts[source_peak_order, , drop = FALSE]
  rownames(matrix) <- unname(source_to_target[source_peak_order])
  list(
    matrix = matrix,
    overlap = overlap,
    genome_qc = list(snp_query = prepared$qc, peak_output = peak_lift$qc)
  )
}

align_genotype_to_snp_matrix <- function(gwas, genotype, snp_matrix, columns) {
  idx_gwas <- match(rownames(snp_matrix), gwas$coord_id)
  if (anyNA(idx_gwas)) stop("SNP accessibility rows do not all match GWAS coordinates.", call. = FALSE)
  gwas_use <- gwas[idx_gwas, , drop = FALSE]
  idx_genotype <- match(as.character(gwas_use[[columns$rsid]]), rownames(genotype))
  if (anyNA(idx_genotype)) stop("SNP accessibility rows do not all match reference genotypes.", call. = FALSE)
  genotype <- genotype[idx_genotype, , drop = FALSE]
  rownames(genotype) <- gwas_use$coord_id
  genotype
}

#' Aggregate library-size-adjusted accessibility by cell group
#'
#' @param feature_matrix Sparse feature-by-cell matrix.
#' @param groups Named group vector aligned to matrix columns.
#' @param cell_scale Optional named per-cell scaling vector. Defaults to one.
#' @return A list with group-specific feature vectors and their unweighted mean baseline.
#' @export
build_group_accessibility <- function(feature_matrix, groups, cell_scale = NULL) {
  require_namespaces("Matrix", "group accessibility aggregation")
  validate_named_matrix(feature_matrix, "feature accessibility")
  if (is.null(names(groups))) names(groups) <- colnames(feature_matrix)
  groups <- groups[colnames(feature_matrix)]
  if (is.null(cell_scale)) cell_scale <- setNames(rep(1, ncol(feature_matrix)), colnames(feature_matrix))
  scale <- as.numeric(cell_scale[colnames(feature_matrix)])
  keep <- !is.na(groups) & is.finite(scale) & scale > 0
  if (!any(keep)) stop("No cells remain after group and library-size alignment.", call. = FALSE)
  matrix <- feature_matrix[, keep, drop = FALSE]
  groups <- as.character(groups[keep])
  group_levels <- sort(unique(groups))
  scaled <- matrix %*% Matrix::Diagonal(x = scale[keep])
  by_group <- lapply(group_levels, function(group) {
    value <- Matrix::rowSums(scaled[, groups == group, drop = FALSE])
    names(value) <- rownames(matrix)
    value
  })
  names(by_group) <- group_levels
  baseline <- rowMeans(do.call(cbind, by_group), na.rm = TRUE)
  names(baseline) <- rownames(matrix)
  list(by_group = by_group, baseline = baseline)
}

compute_cell_scale <- function(atac, assay = "ATAC") {
  require_namespaces(c("Signac", "GenomicRanges", "Matrix"), "library-size scaling")
  counts <- get_atac_counts(atac, assay)
  ranges <- Signac::StringToGRanges(rownames(counts), sep = c("-", "-"))
  autosomal <- chr_to_order(as.character(GenomicRanges::seqnames(ranges))) %in% 1:22
  library_size <- Matrix::colSums(counts[autosomal, , drop = FALSE])
  library_size[library_size <= 0] <- NA_real_
  scale <- stats::median(library_size, na.rm = TRUE) / library_size
  setNames(as.numeric(scale), names(library_size))
}

#' Construct SNP and peak accessibility matrices
#'
#' @param config Parsed EPIC2 configuration.
#' @param traits Optional trait subset.
#' @return Named list of generated files.
#' @export
run_accessibility_stage <- function(config, traits = names(config$gwas$traits)) {
  output_dir <- stage_output_dir(config, "02_accessibility")
  atac_path <- file.path(config$output_dir, "01_atac", "atac.qc.rds")
  atac <- load_serialized(atac_path)
  outputs <- setNames(vector("list", length(traits)), traits)

  for (trait in traits) {
    gwas <- load_serialized(file.path(config$output_dir, "02_gwas", paste0(trait, ".gwas.rds")))
    genotype <- load_serialized(file.path(config$output_dir, "02_gwas", paste0(trait, ".genotype.rds")))
    columns <- deep_merge(config$gwas$columns, trait_config(config, trait)$columns %||% list())
    target_genome <- normalize_genome_build(config$genome$target)
    atac_genome <- normalize_genome_build(config$atac$genome)
    input_builds <- unique(as.character(gwas$EPIC2_input_genome %||% character()))
    input_builds <- input_builds[!is.na(input_builds) & nzchar(input_builds)]
    original_matches_atac <- length(input_builds) == 1L &&
      normalize_genome_build(input_builds) == atac_genome
    target_to_atac_chain <- if (original_matches_atac) NULL else {
      resolve_chain_file(config, target_genome, atac_genome)
    }
    atac_to_target_chain <- resolve_chain_file(config, atac_genome, target_genome)
    snp <- build_snp_accessibility(
      atac, gwas, columns, config$atac$assay,
      workers = config$parameters$workers,
      chunks = config$accessibility$chunks,
      gwas_genome = target_genome,
      atac_genome = atac_genome,
      chain_file = target_to_atac_chain,
      unmapped_action = config$genome$unmapped_action,
      multimapping_action = config$genome$multimapping_action,
      duplicate_action = config$genome$duplicate_action
    )
    peak <- build_peak_accessibility(
      atac, gwas, columns, config$atac$assay,
      gwas_genome = target_genome,
      atac_genome = atac_genome,
      gwas_to_atac_chain = target_to_atac_chain,
      atac_to_gwas_chain = atac_to_target_chain,
      unmapped_action = config$genome$unmapped_action,
      multimapping_action = config$genome$multimapping_action,
      duplicate_action = config$genome$duplicate_action
    )
    aligned_genotype <- align_genotype_to_snp_matrix(gwas, genotype, snp, columns)

    snp_path <- file.path(output_dir, paste0(trait, ".snp_by_cell.rds"))
    peak_path <- file.path(output_dir, paste0(trait, ".peak_by_cell.rds"))
    overlap_path <- file.path(output_dir, paste0(trait, ".snp_peak_overlap.tsv"))
    genome_qc_path <- file.path(output_dir, paste0(trait, ".genome_qc.tsv"))
    genotype_path <- file.path(output_dir, paste0(trait, ".snp_genotype.rds"))
    save_rds_atomic(snp, snp_path)
    save_rds_atomic(peak$matrix, peak_path)
    save_rds_atomic(aligned_genotype, genotype_path)
    write_table(peak$overlap, overlap_path)
    genome_qc <- rbind(
      cbind(step = "snp_query", attr(snp, "genome_qc"), stringsAsFactors = FALSE),
      cbind(step = "peak_snp_query", peak$genome_qc$snp_query, stringsAsFactors = FALSE),
      cbind(step = "peak_output", peak$genome_qc$peak_output, stringsAsFactors = FALSE)
    )
    write_table(genome_qc, genome_qc_path)
    outputs[[trait]] <- list(
      snp = snp_path, peak = peak_path, overlap = overlap_path,
      genotype = genotype_path, genome_qc = genome_qc_path
    )
    message(
      "Accessibility stage complete for ", trait, ": ", nrow(snp),
      " SNPs and ", nrow(peak$matrix), " peaks"
    )
  }
  outputs
}
