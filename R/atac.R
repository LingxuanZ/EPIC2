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
                                  validate = TRUE, barcode_column = NULL) {
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
    if (!is.null(barcode_column) && !barcode_column %in% colnames(metadata)) {
      stop("Fragment barcode metadata column not found: ", barcode_column, call. = FALSE)
    }
    barcodes <- if (is.null(barcode_column)) cells else as.character(metadata[cells, barcode_column])
    if (anyNA(barcodes) || any(!nzchar(barcodes)) || anyDuplicated(barcodes)) {
      stop("Each fragment file requires unique, non-missing raw barcodes.", call. = FALSE)
    }
    Signac::CreateFragmentObject(
      path = path,
      cells = setNames(barcodes, cells),
      validate.fragments = validate,
      verbose = TRUE
    )
  })
  Signac::Fragments(atac[[assay]]) <- fragments
  atac
}

count_fragments_in_peaks <- function(path, cell_map, peak_union, chunk_size = 100000L) {
  path <- assert_file(path, "Fragment file")
  if (is.null(names(cell_map)) || anyNA(names(cell_map)) || any(!nzchar(names(cell_map))) ||
      anyDuplicated(names(cell_map)) || anyNA(cell_map) || any(!nzchar(cell_map)) ||
      anyDuplicated(unname(cell_map))) {
    stop("Fragment barcode mappings must have unique, non-empty cell names and raw barcodes.", call. = FALSE)
  }
  total <- in_peaks <- setNames(numeric(length(cell_map)), names(cell_map))
  if (!length(cell_map)) return(list(total = total, in_peaks = in_peaks))
  peak_index <- split(seq_along(peak_union), as.character(GenomicRanges::seqnames(peak_union)))
  starts <- lapply(peak_index, function(index) GenomicRanges::start(peak_union)[index])
  ends <- lapply(peak_index, function(index) GenomicRanges::end(peak_union)[index])
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  repeat {
    lines <- readLines(connection, n = chunk_size, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    if (!length(lines)) next
    fields <- strsplit(lines, "\t", fixed = TRUE)
    if (any(!lengths(fields) %in% c(5L, 6L))) {
      stop("Fragment records must have five or six tab-separated columns: ", path, call. = FALSE)
    }
    barcodes <- vapply(fields, `[[`, character(1), 4L)
    cell_index <- match(barcodes, unname(cell_map))
    selected <- which(!is.na(cell_index))
    if (!length(selected)) next
    fields <- fields[selected]
    cell_index <- cell_index[selected]
    chromosome <- canonical_chr(vapply(fields, `[[`, character(1), 1L))
    start <- suppressWarnings(as.numeric(vapply(fields, `[[`, character(1), 2L)))
    end <- suppressWarnings(as.numeric(vapply(fields, `[[`, character(1), 3L)))
    if (any(!is.finite(start) | !is.finite(end) | start < 0 | end <= start |
            start != floor(start) | end != floor(end))) {
      stop("Fragment records must use valid zero-based, half-open integer coordinates: ", path,
           call. = FALSE)
    }
    total <- total + tabulate(cell_index, nbins = length(cell_map))
    overlap <- rep(FALSE, length(cell_index))
    for (chr in intersect(unique(chromosome), names(starts))) {
      index <- which(chromosome == chr)
      # BED [start, end) is [start + 1, end] in the one-based peak ranges.
      # Reduced peaks are sorted and disjoint. The last peak starting at or
      # before end overlaps precisely when its inclusive end exceeds start.
      last_peak <- findInterval(end[index], starts[[chr]])
      candidate <- which(last_peak > 0L)
      overlap[index[candidate]] <- ends[[chr]][last_peak[candidate]] > start[index[candidate]]
    }
    in_peaks <- in_peaks + tabulate(cell_index[overlap], nbins = length(cell_map))
  }
  list(total = total, in_peaks = in_peaks)
}

count_atac_fragments <- function(atac, regions, assay) {
  fragments <- Signac::Fragments(atac[[assay]])
  if (!length(fragments)) stop("No fragment objects are attached to the ATAC assay.", call. = FALSE)
  regions <- GenomicRanges::reduce(normalize_granges_seqnames(regions), ignore.strand = TRUE)
  cells <- colnames(atac)
  maps <- list()
  for (fragment in fragments) {
    path <- assert_file(Signac::GetFragmentData(fragment, slot = "path"), "Fragment file")
    cell_map <- Signac::GetFragmentData(fragment, slot = "cells")
    if (is.null(cell_map)) {
      if (length(fragments) > 1L) stop("Multiple fragment objects require explicit cell-to-barcode mappings.", call. = FALSE)
      cell_map <- setNames(cells, cells)
    }
    if (is.null(names(cell_map))) names(cell_map) <- cell_map
    cell_map <- cell_map[names(cell_map) %in% cells]
    maps[[path]] <- c(maps[[path]], cell_map)
  }
  total <- in_peaks <- setNames(numeric(length(cells)), cells)
  for (path in names(maps)) {
    cell_map <- maps[[path]]
    # The same indexed file can be attached more than once after a merge.
    pairs <- unique(data.frame(cell = names(cell_map), barcode = unname(cell_map)))
    cell_map <- setNames(pairs$barcode, pairs$cell)
    if (!length(cell_map)) next
    counted <- count_fragments_in_peaks(path, cell_map, regions)
    mapped <- names(counted$total)
    total[mapped] <- total[mapped] + counted$total
    in_peaks[mapped] <- in_peaks[mapped] + counted$in_peaks
  }
  list(total = total, in_peaks = in_peaks)
}

#' Add fraction-of-fragments-in-peaks metrics
#'
#' Count each fragment-file record once in the denominator and once in the
#' numerator if it overlaps at least one retained peak. The fifth-column read
#' multiplicity is not used. Fragment files are streamed in bounded blocks;
#' barcode mappings are applied separately for each fragment object.
#'
#' @param atac A Seurat object with Signac fragment objects attached.
#' @param assay ATAC assay name.
#' @param genome_sep Two separators used in peak row names.
#' @return The input object with fragment totals and `FRIP` metadata columns.
#'   Cells without fragment records have `NA` FRIP.
#' @export
add_frip <- function(atac, assay = "ATAC", genome_sep = c("-", "-")) {
  require_namespaces(c("Signac", "GenomicRanges"), "FRIP calculation")
  peaks <- Signac::StringToGRanges(rownames(atac[[assay]]), sep = genome_sep)
  counted <- count_atac_fragments(atac, peaks, assay)
  total <- counted$total
  in_peaks <- counted$in_peaks
  ratio <- rep(NA_real_, ncol(atac))
  observed <- total > 0
  ratio[observed] <- in_peaks[observed] / total[observed]
  atac$total_unique_fragments <- as.numeric(total)
  atac$unique_fragments_in_peaks <- as.numeric(in_peaks)
  atac$FRIP <- ratio
  atac
}

#' Add blacklist-overlap ratios
#'
#' Count each fragment record at most once across all blacklist intervals,
#' using the same cell mapping and denominator as [add_frip()].
#' @param atac A Seurat object with fragment objects.
#' @param blacklist A `GRanges` blacklist. Defaults to the Signac blacklist for the explicitly declared genome.
#' @param assay ATAC assay name.
#' @param genome ATAC genome build used to select the default blacklist.
#' @return The input object with blacklist count and ratio metadata columns.
#' @export
add_blacklist_ratio <- function(atac, blacklist = NULL, assay = "ATAC", genome = NULL) {
  require_namespaces(c("Signac", "GenomicRanges"), "blacklist ratio calculation")
  blacklist <- blacklist %||% load_genome_blacklist(genome)
  counted <- count_atac_fragments(atac, blacklist, assay)
  ratio <- rep(NA_real_, ncol(atac))
  observed <- counted$total > 0
  ratio[observed] <- counted$in_peaks[observed] / counted$total[observed]
  atac$total_unique_fragments <- as.numeric(counted$total)
  atac$blacklist_fragments <- as.numeric(counted$in_peaks)
  atac$BlacklistRatio <- ratio
  atac
}

atac_initial_qc <- function(metadata, thresholds) {
  is.finite(metadata$TSS.enrichment) & metadata$TSS.enrichment >= thresholds$tss_min &
    is.finite(metadata$FRIP) & metadata$FRIP >= thresholds$frip_min & metadata$FRIP <= 1
}

atac_count_cutoff <- function(metadata, count_quantile_min, assay = "ATAC",
                              thresholds = EPIC2_defaults()$atac$qc) {
  count_column <- paste0("nCount_", assay)
  if (!count_column %in% colnames(metadata)) {
    stop("Missing ATAC metadata required for cell QC: ", count_column, call. = FALSE)
  }
  keep <- atac_initial_qc(metadata, thresholds) & is.finite(metadata[[count_column]])
  if (!any(keep)) stop("No cells pass TSS/FRIP QC with finite counts.", call. = FALSE)
  stats::quantile(
    metadata[[count_column]][keep],
    probs = count_quantile_min,
    na.rm = TRUE,
    names = FALSE
  )
}

#' Filter nuclei using configured QC thresholds
#'
#' @param atac A Seurat/Signac object.
#' @param thresholds Named list of cell-level thresholds.
#' @param group_column Metadata column used for minimum group size filtering.
#' @param assay Chromatin assay name; cell counts are read from `nCount_<assay>`.
#' @return The filtered object.
#' @export
filter_atac_cells <- function(atac, thresholds = list(), group_column = "tissue_cell_type",
                              assay = "ATAC") {
  defaults <- EPIC2_defaults()$atac$qc
  thresholds <- deep_merge(defaults, thresholds)
  metadata <- get_atac_metadata(atac)
  count_column <- paste0("nCount_", assay)
  required <- c("TSS.enrichment", "FRIP", count_column, "BlacklistRatio", "nucleosome_signal", group_column)
  missing <- setdiff(required, colnames(metadata))
  if (length(missing)) {
    stop("Missing ATAC metadata required for cell QC: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  count_cutoff <- atac_count_cutoff(metadata, thresholds$count_quantile_min, assay, thresholds)
  counts <- metadata[[count_column]]
  keep <- atac_initial_qc(metadata, thresholds) &
    is.finite(counts) & counts >= count_cutoff &
    is.finite(metadata$BlacklistRatio) & metadata$BlacklistRatio >= 0 & metadata$BlacklistRatio < thresholds$blacklist_max &
    is.finite(metadata$nucleosome_signal) & metadata$nucleosome_signal <= thresholds$nucleosome_max

  groups <- as.character(metadata[[group_column]])
  group_sizes <- table(groups[keep & !is.na(groups)])
  allowed <- names(group_sizes)[group_sizes >= thresholds$min_cells_per_group]
  keep <- keep & !is.na(groups) & groups %in% allowed
  if (!any(keep)) stop("Cell QC removed every cell; inspect configured thresholds.", call. = FALSE)
  atac[, rownames(metadata)[keep]]
}

#' Filter ATAC peaks using configured QC thresholds
#'
#' @param atac A Seurat/Signac object.
#' @param thresholds Named list containing peak width thresholds.
#' @param assay ATAC assay name.
#' @param blacklist A `GRanges` blacklist. Defaults to the Signac blacklist for the explicitly declared genome.
#' @param genome ATAC genome build used to select the default blacklist.
#' @return The filtered object.
#' @export
filter_atac_peaks <- function(atac, thresholds = list(), assay = "ATAC", blacklist = NULL,
                              genome = NULL) {
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
  config <- normalize_common_genome_config(config)
  lock <- acquire_stage_lock(config, "atac")
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  output_dir <- stage_output_dir(config, "01_atac")
  input <- resolve_config_path(config, config$atac$input)
  atac <- load_serialized(input, config$atac$object_name %||% NULL)
  assay <- config$atac$assay
  if (!assay %in% names(atac)) stop("ATAC assay not found in input object: ", assay, call. = FALSE)
  require_namespaces(c("Seurat", "Signac", "GenomicRanges"), "ATAC assay selection")
  Seurat::DefaultAssay(atac) <- assay
  atac_genome <- config$genome$build
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
      validate = config$atac$validate_fragments %||% TRUE,
      barcode_column = config$atac$fragment_barcode_column
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
  count_cutoff <- atac_count_cutoff(metadata_before, config$atac$qc$count_quantile_min,
                                   assay, config$atac$qc)
  before <- data.frame(
    stage = "before", cells = ncol(atac), peaks = nrow(get_atac_counts(atac, assay))
  )
  atac <- filter_atac_cells(atac, config$atac$qc, config$atac$group_column, assay)
  atac <- filter_atac_peaks(atac, config$atac$qc, assay, genome = atac_genome)
  if (isTRUE(config$atac$recompute_embedding)) atac <- recompute_atac_reductions(atac, assay)
  after <- data.frame(
    stage = "after", cells = ncol(atac), peaks = nrow(get_atac_counts(atac, assay))
  )

  output <- file.path(output_dir, "atac.qc.rds")
  save_rds_atomic(atac, output)
  write_table(rbind(before, after), file.path(output_dir, "qc_dimensions.tsv"))
  genome_metadata_qc$genome_build <- atac_genome
  genome_metadata_qc$operation <- "validation_only"
  genome_metadata_qc$conversion_performed <- FALSE
  genome_metadata_qc$coordinate_output <- "Input coordinates retained; all inputs must already use genome.build"
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
    parameter = paste0("derived.nCount_", assay, "_min"),
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
