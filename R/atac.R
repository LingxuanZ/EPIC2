load_genome_blacklist <- function(genome) {
  genome <- normalize_genome_build(genome)
  require_namespaces("Signac", paste0(genome, " blacklist filtering"))
  dataset <- if (genome == "hg19") "blacklist_hg19" else "blacklist_hg38"
  env <- new.env(parent = emptyenv())
  suppressWarnings(utils::data(list = dataset, package = "Signac", envir = env))
  if (!exists(dataset, envir = env, inherits = FALSE)) {
    stop("Signac::", dataset, " is unavailable.", call. = FALSE)
  }
  get(dataset, envir = env, inherits = FALSE)
}

load_hg19_blacklist <- function() {
  load_genome_blacklist("hg19")
}

attach_fragment_files <- function(atac, fragment_files, tissue_column, assay = "ATAC",
                                  validate = TRUE) {
  require_namespaces(c("Signac", "Seurat"), "fragment attachment")
  metadata <- get_atac_metadata(atac)
  if (!tissue_column %in% colnames(metadata)) {
    stop("ATAC metadata column not found: ", tissue_column, call. = FALSE)
  }
  if (is.null(names(fragment_files)) || any(!nzchar(names(fragment_files)))) {
    stop("atac.fragment_files must be a named tissue-to-file mapping.", call. = FALSE)
  }

  fragments <- lapply(names(fragment_files), function(tissue) {
    path <- assert_file(fragment_files[[tissue]], paste0("Fragment file for ", tissue))
    cells <- rownames(metadata)[as.character(metadata[[tissue_column]]) == tissue]
    if (length(cells) == 0L) stop("No cells found for fragment tissue: ", tissue, call. = FALSE)
    Signac::CreateFragmentObject(
      path = path,
      cells = cells,
      validate.fragments = validate,
      verbose = TRUE
    )
  })
  Signac::Fragments(atac[[assay]]) <- fragments
  atac
}

#' Add fraction-of-reads-in-peaks metrics
#'
#' @param atac A Seurat object with Signac fragment objects attached.
#' @param assay ATAC assay name.
#' @param genome_sep Two separators used in peak row names.
#' @return The input object with fragment totals and `FRIP` metadata columns.
#' @export
add_frip <- function(atac, assay = "ATAC", genome_sep = c("-", "-")) {
  require_namespaces(c("Signac", "GenomicRanges", "Matrix"), "FRIP calculation")
  fragments <- Signac::Fragments(atac[[assay]])
  if (length(fragments) == 0L) stop("No fragment objects are attached to the ATAC assay.", call. = FALSE)

  totals <- lapply(fragments, function(fragment) {
    Signac::CountFragments(fragment@path, cells = names(fragment@cells))
  })
  totals <- do.call(rbind, totals)
  total_by_cell <- rowsum(totals$frequency_count, group = totals$CB, reorder = FALSE)[, 1L]

  peaks <- Signac::StringToGRanges(rownames(atac[[assay]]), sep = genome_sep)
  peak_union <- GenomicRanges::reduce(peaks, ignore.strand = TRUE)
  peak_matrix <- Signac::FeatureMatrix(
    fragments = fragments,
    features = peak_union,
    cells = colnames(atac)
  )
  in_peaks <- Matrix::colSums(peak_matrix)
  total_aligned <- total_by_cell[colnames(atac)]
  atac$total_unique_fragments <- as.numeric(total_aligned)
  atac$unique_fragments_in_peaks <- as.numeric(in_peaks[colnames(atac)])
  atac$FRIP <- atac$unique_fragments_in_peaks / atac$total_unique_fragments
  atac
}

#' Add blacklist-overlap ratios
#'
#' @param atac A Seurat object with fragment objects and total fragment metadata.
#' @param blacklist A `GRanges` blacklist. Defaults to Signac's hg19 blacklist.
#' @param assay ATAC assay name.
#' @param genome ATAC genome build used to select the default blacklist.
#' @return The input object with blacklist count and ratio metadata columns.
#' @export
add_blacklist_ratio <- function(atac, blacklist = NULL, assay = "ATAC", genome = "hg19") {
  require_namespaces(c("Signac", "GenomicRanges", "Matrix"), "blacklist ratio calculation")
  metadata <- get_atac_metadata(atac)
  if (!"total_unique_fragments" %in% colnames(metadata)) {
    stop("Run add_frip() before add_blacklist_ratio().", call. = FALSE)
  }
  blacklist <- blacklist %||% load_genome_blacklist(genome)
  blacklist <- GenomicRanges::reduce(blacklist, ignore.strand = TRUE)
  matrix <- Signac::FeatureMatrix(
    fragments = Signac::Fragments(atac[[assay]]),
    features = blacklist,
    cells = colnames(atac)
  )
  counts <- Matrix::colSums(matrix)
  atac$blacklist_fragments <- as.numeric(counts[colnames(atac)])
  atac$BlacklistRatio <- atac$blacklist_fragments / atac$total_unique_fragments
  atac
}

atac_count_cutoff <- function(metadata, count_quantile_min) {
  stats::quantile(
    metadata$nCount_ATAC,
    probs = count_quantile_min,
    na.rm = TRUE,
    names = FALSE
  )
}

#' Filter nuclei using manuscript QC thresholds
#'
#' @param atac A Seurat/Signac object.
#' @param thresholds Named list of cell-level thresholds.
#' @param group_column Metadata column used for minimum group size filtering.
#' @return The filtered object.
#' @export
filter_atac_cells <- function(atac, thresholds = list(), group_column = "tissue_cell_type") {
  defaults <- EPIC2_defaults()$atac$qc
  thresholds <- deep_merge(defaults, thresholds)
  metadata <- get_atac_metadata(atac)
  required <- c("TSS.enrichment", "FRIP", "nCount_ATAC", "BlacklistRatio", "nucleosome_signal", group_column)
  missing <- setdiff(required, colnames(metadata))
  if (length(missing)) {
    stop("Missing ATAC metadata required for cell QC: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  count_cutoff <- atac_count_cutoff(metadata, thresholds$count_quantile_min)
  keep <- is.finite(metadata$TSS.enrichment) & metadata$TSS.enrichment >= thresholds$tss_min &
    is.finite(metadata$FRIP) & metadata$FRIP >= thresholds$frip_min &
    is.finite(metadata$nCount_ATAC) & metadata$nCount_ATAC >= count_cutoff &
    is.finite(metadata$BlacklistRatio) & metadata$BlacklistRatio < thresholds$blacklist_max &
    is.finite(metadata$nucleosome_signal) & metadata$nucleosome_signal <= thresholds$nucleosome_max

  groups <- as.character(metadata[[group_column]])
  group_sizes <- table(groups[keep & !is.na(groups)])
  allowed <- names(group_sizes)[group_sizes >= thresholds$min_cells_per_group]
  keep <- keep & !is.na(groups) & groups %in% allowed
  if (!any(keep)) stop("Cell QC removed every cell; inspect configured thresholds.", call. = FALSE)
  atac[, rownames(metadata)[keep]]
}

#' Filter ATAC peaks using manuscript QC thresholds
#'
#' @param atac A Seurat/Signac object.
#' @param thresholds Named list containing peak width thresholds.
#' @param assay ATAC assay name.
#' @param blacklist A `GRanges` blacklist. Defaults to Signac's hg19 blacklist.
#' @param genome ATAC genome build used to select the default blacklist.
#' @return The filtered object.
#' @export
filter_atac_peaks <- function(atac, thresholds = list(), assay = "ATAC", blacklist = NULL,
                              genome = "hg19") {
  require_namespaces(c("Signac", "GenomicRanges", "Matrix"), "peak QC")
  thresholds <- deep_merge(EPIC2_defaults()$atac$qc, thresholds)
  counts <- get_atac_counts(atac, assay)
  peaks <- Signac::StringToGRanges(rownames(counts), sep = c("-", "-"))
  peaks <- normalize_granges_seqnames(peaks)

  keep <- Matrix::rowSums(counts) >= thresholds$peak_min_total_counts &
    Matrix::rowSums(counts > 0) >= thresholds$peak_min_cells &
    GenomicRanges::width(peaks) >= thresholds$peak_width_min &
    GenomicRanges::width(peaks) <= thresholds$peak_width_max
  if (isTRUE(thresholds$peak_autosomes_only)) {
    keep <- keep & chr_to_order(as.character(GenomicRanges::seqnames(peaks))) %in% 1:22
  }
  if (isTRUE(thresholds$peak_remove_blacklist)) {
    blacklist <- blacklist %||% load_genome_blacklist(genome)
    keep <- keep & !IRanges::overlapsAny(peaks, blacklist, ignore.strand = TRUE)
  }
  if (!any(keep)) stop("Peak QC removed every peak; inspect coordinates and thresholds.", call. = FALSE)
  atac[rownames(counts)[keep], ]
}

recompute_atac_reductions <- function(atac, assay = "ATAC") {
  require_namespaces(c("Signac", "Seurat"), "ATAC dimensionality reduction")
  Seurat::DefaultAssay(atac) <- assay
  atac <- Signac::RunTFIDF(atac)
  atac <- Signac::FindTopFeatures(atac, min.cutoff = "q0")
  atac <- Signac::RunSVD(atac)
  atac <- Seurat::RunUMAP(atac, reduction = "lsi", dims = 2:30)
  atac <- Seurat::FindNeighbors(atac, reduction = "lsi", dims = 2:30)
  Seurat::FindClusters(atac, verbose = FALSE, algorithm = 3)
}

#' ATAC preparation and QC
#'
#' @param config Parsed EPIC2 configuration.
#' @return Path to the filtered ATAC object.
#' @export
run_atac_stage <- function(config) {
  output_dir <- stage_output_dir(config, "01_atac")
  input <- resolve_config_path(config, config$atac$input)
  atac <- load_serialized(input, config$atac$object_name %||% NULL)
  assay <- config$atac$assay
  if (!assay %in% names(atac)) stop("ATAC assay not found in input object: ", assay, call. = FALSE)
  require_namespaces(c("Seurat", "Signac", "GenomicRanges"), "ATAC assay selection")
  Seurat::DefaultAssay(atac) <- assay
  atac_genome <- normalize_genome_build(config$atac$genome)
  target_genome <- normalize_genome_build(config$genome$target)
  input_ranges <- tryCatch(
    GenomicRanges::granges(atac[[assay]]),
    error = function(e) Signac::StringToGRanges(rownames(get_atac_counts(atac, assay)), sep = c("-", "-"))
  )
  genome_metadata_qc <- check_range_genome_declaration(input_ranges, atac_genome, "ATAC assay")
  validate_genome_coordinates(
    as.character(GenomicRanges::seqnames(input_ranges)), GenomicRanges::start(input_ranges),
    GenomicRanges::end(input_ranges), atac_genome, "ATAC assay peaks"
  )

  if (isTRUE(config$atac$attach_fragments)) {
    files <- config$atac$fragment_files
    files <- setNames(
      vapply(files, function(path) resolve_config_path(config, path), character(1)),
      names(files)
    )
    atac <- attach_fragment_files(
      atac, files, config$atac$tissue_column, assay,
      validate = config$atac$validate_fragments %||% TRUE
    )
  }

  if (isTRUE(config$atac$compute_qc_metrics)) {
    require_namespaces("Signac", "ATAC QC metrics")
    atac <- Signac::NucleosomeSignal(atac)
    atac <- Signac::TSSEnrichment(atac)
    atac <- add_frip(atac, assay)
    atac <- add_blacklist_ratio(atac, assay = assay, genome = atac_genome)
  }

  metadata_before <- get_atac_metadata(atac)
  count_cutoff <- atac_count_cutoff(metadata_before, config$atac$qc$count_quantile_min)
  before <- data.frame(
    stage = "before", cells = ncol(atac), peaks = nrow(get_atac_counts(atac, assay))
  )
  atac <- filter_atac_cells(atac, config$atac$qc, config$atac$group_column)
  atac <- filter_atac_peaks(atac, config$atac$qc, assay, genome = atac_genome)
  if (isTRUE(config$atac$recompute_embedding)) atac <- recompute_atac_reductions(atac, assay)
  after <- data.frame(
    stage = "after", cells = ncol(atac), peaks = nrow(get_atac_counts(atac, assay))
  )

  output <- file.path(output_dir, "atac.qc.rds")
  save_rds_atomic(atac, output)
  write_table(rbind(before, after), file.path(output_dir, "qc_dimensions.tsv"))
  genome_metadata_qc$target_genome <- target_genome
  genome_metadata_qc$coordinate_output <- if (atac_genome == target_genome) {
    "ATAC and downstream outputs use the same build"
  } else {
    "ATAC fragments remain in the input build; downstream SNP/peak IDs are lifted"
  }
  write_table(genome_metadata_qc, file.path(output_dir, "genome_qc.tsv"))
  qc_values <- vapply(config$atac$qc, function(value) {
    if (length(value) == 0L) "null" else paste(as.character(value), collapse = ",")
  }, character(1))
  qc_manifest <- data.frame(
    parameter = paste0("atac.qc.", names(qc_values)),
    value = unname(qc_values),
    stringsAsFactors = FALSE
  )
  qc_manifest <- rbind(qc_manifest, data.frame(
    parameter = "derived.nCount_ATAC_min",
    value = format(count_cutoff, scientific = FALSE, trim = TRUE),
    stringsAsFactors = FALSE
  ))
  write_table(qc_manifest, file.path(output_dir, "qc_thresholds.tsv"))

  metadata <- get_atac_metadata(atac)
  group_summary <- as.data.frame(table(metadata[[config$atac$group_column]], useNA = "ifany"))
  colnames(group_summary) <- c("group", "n_cells")
  write_table(group_summary, file.path(output_dir, "group_cell_counts.tsv"))
  message("ATAC stage complete: ", output)
  output
}
