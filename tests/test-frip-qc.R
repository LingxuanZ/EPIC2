# Exercise FRIP against real indexed BED fragments and valid Signac mappings.
source_directory <- Sys.getenv("EPIC2_MODIFIED_R_DIR", "R")
if (dir.exists(source_directory)) {
  epic2 <- new.env(parent = globalenv())
  for (path in sort(list.files(source_directory, pattern = "\\.R$", full.names = TRUE))) {
    sys.source(path, envir = epic2)
  }
} else {
  epic2 <- asNamespace("EPIC2")
}
required <- c("Signac", "Seurat", "Rsamtools", "Matrix", "GenomicRanges", "IRanges", "GenomeInfoDb", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("FRIP regression test requires: ", paste(missing, collapse = ", "))

run_frip_regression <- function() {
  directory <- tempfile("EPIC2-frip-qc-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  peaks <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100L, 500L), c(200L, 600L)))
  GenomeInfoDb::genome(peaks) <- "hg19"
  peak_ids <- paste0("chr1-", GenomicRanges::start(peaks), "-", GenomicRanges::end(peaks))
  make_fragment <- function(name, records, mapping) {
    path <- file.path(directory, paste0(name, ".tsv"))
    records <- records[order(records$chr, records$start, records$end), , drop = FALSE]
    writeLines("# BED coordinates; fifth column is read multiplicity", path)
    utils::write.table(records, path, sep = "\t", quote = FALSE, row.names = FALSE,
                       col.names = FALSE, append = TRUE)
    compressed <- Rsamtools::bgzip(path, overwrite = TRUE)
    Rsamtools::indexTabix(compressed, format = "bed")
    Signac::CreateFragmentObject(compressed, cells = mapping,
                                  validate.fragments = TRUE, verbose = FALSE)
  }
  first <- data.frame(
    chr = "chr1",
    # For AAAC-1: no overlap, start boundary, both peaks, end boundary,
    # just after first peak, just after last peak. Three of six overlap.
    start = c(89L, 99L, 150L, 199L, 200L, 600L, 110L, 110L),
    end = c(99L, 100L, 550L, 200L, 201L, 601L, 120L, 120L),
    cell = c(rep("AAAC-1", 6L), "AAAG-1", "unselected"),
    count = c(1L, 7L, 50L, 2L, 1L, 1L, 9L, 1L)
  )
  second <- data.frame(
    chr = "chr1", start = c(250L, 510L, 700L), end = c(260L, 520L, 710L),
    cell = c("AAAC-1", "AAAC-1", "AAAG-1"), count = c(3L, 12L, 4L)
  )
  maps <- list(c(sample1_AAAC = "AAAC-1", sample1_AAAG = "AAAG-1"),
               c(sample2_AAAC = "AAAC-1", sample2_AAAG = "AAAG-1"))
  fragments <- list(make_fragment("sample1", first, maps[[1L]]),
                    make_fragment("sample2", second, maps[[2L]]))
  cells <- c(names(maps[[1L]]), names(maps[[2L]]))
  make_atac <- function(assay, fragments, cells) {
    counts <- matrix(rep(seq_along(cells) * 5, each = 2L), nrow = 2L,
                      dimnames = list(peak_ids, cells))
    chromatin <- Signac::CreateChromatinAssay(
      counts = Matrix::Matrix(counts, sparse = TRUE), ranges = peaks, fragments = fragments,
      min.cells = 0, min.features = 0, validate.fragments = FALSE
    )
    Seurat::CreateSeuratObject(counts = chromatin, assay = assay)
  }
  atac <- make_atac("ATAC", fragments, cells)
  actual <- epic2$add_frip(atac)
  stopifnot(
    identical(as.numeric(actual$total_unique_fragments), c(6, 1, 2, 1)),
    identical(as.numeric(actual$unique_fragments_in_peaks), c(3, 1, 1, 0)),
    identical(as.numeric(actual$FRIP), c(0.5, 1, 0.5, 0)),
    all(actual$FRIP >= 0 & actual$FRIP <= 1)
  )
  blacklisted <- epic2$add_blacklist_ratio(atac, blacklist = peaks)
  stopifnot(identical(as.numeric(blacklisted$BlacklistRatio), c(0.5, 1, 0.5, 0)),
            all(blacklisted$BlacklistRatio >= 0 & blacklisted$BlacklistRatio <= 1))
  # Force boundaries between input blocks, including a comment-only block.
  chunked <- epic2$count_fragments_in_peaks(
    Signac::GetFragmentData(fragments[[1L]], "path"), maps[[1L]],
    GenomicRanges::reduce(peaks, ignore.strand = TRUE), chunk_size = 1L
  )
  stopifnot(identical(unname(chunked$total), c(6, 1)),
            identical(unname(chunked$in_peaks), c(3, 1)))
  # MT and chrM are the same standard mitochondrial chromosome label.
  mitochondrial_path <- file.path(directory, "mitochondrial.tsv")
  writeLines("MT\t99\t120\tAAAC-1\t1", mitochondrial_path)
  mitochondrial <- epic2$count_fragments_in_peaks(mitochondrial_path,
    c(cell1 = "AAAC-1"), GenomicRanges::GRanges("chrM", IRanges::IRanges(100, 200)))
  stopifnot(mitochondrial$total[[1]] == 1, mitochondrial$in_peaks[[1]] == 1)

  # Ordinary barcodes with no renaming retain the same fragment-based definition.
  standard_records <- data.frame(
    chr = c("chr1", "chr1", "chr2"), start = c(110L, 510L, 110L),
    end = c(120L, 520L, 120L), cell = c("cell1", "cell2", "cell1"), count = c(8L, 2L, 1L)
  )
  standard_cells <- c("cell1", "cell2")
  standard_fragment <- make_fragment("standard", standard_records,
                                     setNames(standard_cells, standard_cells))
  standard <- epic2$add_frip(make_atac("ATAC", list(standard_fragment), standard_cells))
  stopifnot(identical(as.numeric(standard$total_unique_fragments), c(2, 1)),
            identical(as.numeric(standard$unique_fragments_in_peaks), c(1, 1)),
            identical(as.numeric(standard$FRIP), c(0.5, 1)))

  # Cells absent from the mapping or with a mapped but unobserved barcode
  # have zero counts and NA FRIP, rather than an infinite ratio.
  missing_cell <- "unobserved"
  no_mapping <- epic2$add_frip(make_atac(
    "ATAC", list(standard_fragment), c(standard_cells, missing_cell)
  ))
  absent_fragment <- Signac::CreateFragmentObject(
    Signac::GetFragmentData(standard_fragment, "path"),
    cells = setNames(c(standard_cells, missing_cell), c(standard_cells, missing_cell)),
    validate.fragments = FALSE, verbose = FALSE
  )
  no_records <- epic2$add_frip(make_atac(
    "ATAC", list(absent_fragment), c(standard_cells, missing_cell)
  ))
  for (zero_case in list(no_mapping, no_records)) {
    stopifnot(identical(as.numeric(zero_case$total_unique_fragments), c(2, 1, 0)),
              identical(as.numeric(zero_case$unique_fragments_in_peaks), c(1, 1, 0)),
              is.na(zero_case$FRIP[[3L]]), !is.infinite(zero_case$FRIP[[3L]]))
  }

  # Custom assay QC must use nCount_peaks even if stale nCount_ATAC also exists.
  custom <- make_atac("peaks", fragments, cells)
  custom$TSS.enrichment <- c(5, 6, 7, 8)
  custom$FRIP <- 0.8
  custom$BlacklistRatio <- 0
  custom$nucleosome_signal <- 1
  custom$tissue_cell_type <- "group1"
  custom$nCount_ATAC <- c(400, 300, 200, 100)
  thresholds <- list(count_quantile_min = 0.5, min_cells_per_group = 1L,
                     peak_remove_blacklist = FALSE)
  filtered <- epic2$filter_atac_cells(custom, thresholds, assay = "peaks")
  stopifnot(identical(colnames(filtered), cells[3:4]))
  # Absence of an ATAC-specific column must also work.
  custom$nCount_ATAC <- NULL
  filtered <- epic2$filter_atac_cells(custom, thresholds, assay = "peaks")
  stopifnot(identical(colnames(filtered), cells[3:4]))
  settings <- epic2$EPIC2_defaults()$plots
  histogram <- epic2$plot_atac_qc_metrics(custom, settings, assay = "peaks")
  density <- epic2$plot_atac_density_scatter(custom, settings, assay = "peaks")
  count_label <- "log10(nCount_peaks + 1)"
  stopifnot(count_label %in% levels(histogram$data$metric),
            isTRUE(all.equal(histogram$data$value[histogram$data$metric == count_label],
                             log10(custom$nCount_peaks + 1), check.attributes = FALSE)),
            identical(density$labels$x, count_label),
            isTRUE(all.equal(density$data$count, log10(custom$nCount_peaks + 1),
                             check.attributes = FALSE)))
  invisible(ggplot2::ggplot_build(histogram))
  invisible(ggplot2::ggplot_build(density))
  default_histogram <- epic2$plot_atac_qc_metrics(standard, settings)
  stopifnot("log10(nCount_ATAC + 1)" %in% levels(default_histogram$data$metric))

  # The pipeline must propagate the assay to its cutoff, filter, and QC manifest.
  config <- epic2$EPIC2_defaults()
  config$genome$build <- "hg19"
  config$output_dir <- file.path(directory, "output")
  config$atac$input <- file.path(directory, "custom.rds")
  config$atac$assay <- "peaks"
  config$atac$qc <- epic2$deep_merge(config$atac$qc, thresholds)
  saveRDS(custom, config$atac$input)
  result <- epic2$run_atac_stage(config)
  stopifnot(identical(colnames(readRDS(result)), cells[3:4]))
  manifest <- utils::read.delim(file.path(config$output_dir, "01_atac", "qc_thresholds.tsv"))
  stopifnot("derived.nCount_peaks_min" %in% manifest$parameter,
            !"derived.nCount_ATAC_min" %in% manifest$parameter)
}
run_frip_regression()
message("FRIP barcode, overlap, BED boundary, streaming, and custom-assay QC regressions passed.")
