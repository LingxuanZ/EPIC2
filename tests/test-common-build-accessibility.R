# Optional historical comparison: set EPIC2_BASE_R_DIR to an earlier source R/.
# No user data or external reference panels are read by this synthetic test.
load_source_environment <- function(directory) {
  result <- new.env(parent = globalenv())
  for (path in sort(list.files(directory, pattern = "\\.R$", full.names = TRUE))) {
    sys.source(path, envir = result)
  }
  result
}
modified_directory <- Sys.getenv("EPIC2_MODIFIED_R_DIR", file.path("R"))
if (dir.exists(modified_directory)) {
  modified <- load_source_environment(modified_directory)
} else {
  # R CMD check runs tests outside the source tree; use the installed code.
  namespace <- asNamespace("EPIC2")
  modified <- new.env(parent = namespace)
  for (name in ls(namespace, all.names = TRUE)) {
    value <- get(name, envir = namespace, inherits = FALSE)
    if (is.function(value) && identical(environment(value), namespace)) {
      environment(value) <- modified
      assign(name, value, envir = modified)
    }
  }
}
base_directory <- Sys.getenv("EPIC2_BASE_R_DIR", "")
if (nzchar(base_directory) && !dir.exists(base_directory)) {
  stop("EPIC2_BASE_R_DIR must name an existing historical source R directory.")
}
baseline <- if (nzchar(base_directory)) load_source_environment(base_directory) else NULL
if (is.null(baseline)) {
  message("Historical comparison skipped: set EPIC2_BASE_R_DIR to an earlier source R directory.")
}
required <- c("Signac", "Seurat", "Rsamtools", "Matrix", "GenomicRanges", "IRanges",
              "GenomeInfoDb", "S4Vectors", "BiocParallel")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Synthetic accessibility test requires: ", paste(missing, collapse = ", "))

expect_error_message <- function(expression, pattern) {
  message <- tryCatch({ force(expression); NA_character_ }, error = conditionMessage)
  stopifnot(!is.na(message), grepl(pattern, message, ignore.case = TRUE))
}
same_matrix <- function(first, second) {
  stopifnot(
    identical(dim(first), dim(second)), identical(dimnames(first), dimnames(second)),
    isTRUE(all.equal(as.numeric(first), as.numeric(second), tolerance = 1e-12))
  )
}
same_numeric <- function(first, second) {
  stopifnot(
    identical(names(first), names(second)), identical(dim(first), dim(second)),
    identical(dimnames(first), dimnames(second)),
    isTRUE(all.equal(as.numeric(first), as.numeric(second), tolerance = 1e-12))
  )
}

# These helpers must never be called by the revised common-build workflow.
for (name in c("lift_genomic_ranges", "resolve_chain_file", "harmonize_gwas_coordinates")) {
  assign(name, function(...) stop("Unexpected coordinate conversion helper call"), envir = modified)
}

run_synthetic_accessibility_test <- function(build) {
  directory <- tempfile("EPIC2-common-build-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  cells <- paste0("cell", seq_len(8L))
  good_position <- c(1000000L, 1000050L, 1001000L, 1002000L, 1003000L)
  position <- c(good_position, 1001500L, 1001500L, 1003500L)
  peaks <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(
      start = c(999900L, 1000800L, 1001800L, 1002800L, 1001400L, 1003400L),
      end = c(1000100L, 1001200L, 1002200L, 1003200L, 1001600L, 1003600L)
    )
  )
  GenomeInfoDb::genome(peaks) <- build
  fragment_rows <- do.call(rbind, lapply(seq_along(unique(position)), function(index) {
    data.frame(
      chr = "chr1", start = unique(position)[index] - 5L,
      end = unique(position)[index] + 5L, cell = cells,
      count = as.integer(1L + (index + seq_along(cells)) %% 3L)
    )
  }))
  fragment_rows <- fragment_rows[order(fragment_rows$start, fragment_rows$end, fragment_rows$cell), ]
  fragment_path <- file.path(directory, "fragments.tsv")
  utils::write.table(fragment_rows, fragment_path, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  fragment_gz <- Rsamtools::bgzip(fragment_path, overwrite = TRUE)
  Rsamtools::indexTabix(fragment_gz, format = "bed")
  fragments <- Signac::CreateFragmentObject(
    fragment_gz, cells = setNames(cells, cells), validate.fragments = TRUE, verbose = FALSE
  )
  counts <- Signac::FeatureMatrix(fragments, features = peaks, cells = cells, verbose = FALSE)
  chromatin <- Signac::CreateChromatinAssay(
    counts, ranges = peaks, fragments = fragments, min.cells = 0, min.features = 0,
    validate.fragments = FALSE
  )
  metadata <- data.frame(
    tissue = rep(c("A", "B"), each = 4L),
    tissue_cell_type = rep(c("A_type", "B_type"), each = 4L), row.names = cells
  )
  atac <- Seurat::CreateSeuratObject(counts = chromatin, assay = "ATAC", meta.data = metadata)
  gwas <- data.frame(
    rsid = paste0("rs", seq_along(position)), chr = 1L, pos = position,
    A1 = "A", A2 = "G", MAF = 0.25, Beta = seq_along(position) / 10,
    SE = 0.1, P = 0.01, build = build, stringsAsFactors = FALSE
  )
  gwas <- gwas[c(4L, 8L, 1L, 6L, 3L, 7L, 5L, 2L), ]
  genotype <- outer(seq_len(7L), seq_len(12L), function(variant, sample) {
    (variant + sample + variant * as.integer(sample %% 2L == 0L)) %% 3L
  })
  rownames(genotype) <- paste0("rs", seq_len(7L))
  colnames(genotype) <- paste0("reference", seq_len(12L))
  genotype <- genotype[c(7L, 2L, 5L, 1L, 6L, 4L, 3L), , drop = FALSE]
  columns <- list(rsid = "rsid", chr = "chr", pos = "pos", a1 = "A1", a2 = "A2",
                  maf = "MAF", p = "P", beta = "Beta", se = "SE")
  arguments <- list(
    gwas = gwas, genotype = genotype, columns = columns, maf_min = 0.01,
    remove_ambiguous = FALSE, ld_prune = FALSE, genotype_format = "dosage",
    genotype_orientation = "variants_by_samples", genome_column = "build"
  )
  current <- do.call(modified$harmonize_gwas, c(arguments, list(genome_build = build)))
  expected_ids <- paste0("chr1-", good_position, "-", good_position)
  stopifnot(
    identical(current$gwas$rsid, paste0("rs", seq_len(5L))),
    identical(current$gwas$coord_id, expected_ids),
    identical(rownames(current$genotype), current$gwas$rsid),
    current$summary$n_reference_matched == 7L,
    current$summary$n_duplicate_coordinates == 2L,
    current$summary$n_final == 5L
  )
  snp <- modified$build_snp_accessibility(
    atac, current$gwas, columns, workers = 1L, chunks = 2L, genome_build = build
  )
  peak <- modified$build_peak_accessibility(atac, current$gwas, columns, genome_build = build)
  aligned <- modified$align_genotype_to_snp_matrix(current$gwas, current$genotype, snp, columns)
  stopifnot(
    identical(rownames(snp), expected_ids), identical(rownames(aligned), expected_ids),
    nrow(peak$matrix) == 4L, nrow(peak$overlap) == 5L,
    all(peak$overlap$snp %in% rownames(snp)),
    all(peak$overlap$peak %in% rownames(peak$matrix))
  )
  same_matrix(aligned, {
    expected <- genotype[match(current$gwas$rsid, rownames(genotype)), , drop = FALSE]
    rownames(expected) <- expected_ids
    expected
  })
  snp_ld <- modified$compute_snp_ld_covariance(aligned, current$gwas, columns)
  peak_ld <- modified$build_peak_ld_model(current$gwas, aligned, peak$matrix, peak$overlap)
  stopifnot(all(is.finite(as.numeric(snp_ld))), all(is.finite(peak_ld$response)))

  if (!is.null(baseline)) {
    previous <- do.call(baseline$harmonize_gwas, c(arguments, list(source_genome = build)))
    for (column in names(previous$gwas)) {
      stopifnot(identical(previous$gwas[[column]], current$gwas[[column]]))
    }
    same_matrix(previous$genotype, current$genotype)
    previous_snp <- baseline$build_snp_accessibility(
      atac, previous$gwas, columns, workers = 1L, chunks = 2L, gwas_genome = build
    )
    previous_peak <- baseline$build_peak_accessibility(atac, previous$gwas, columns, gwas_genome = build)
    previous_aligned <- baseline$align_genotype_to_snp_matrix(
      previous$gwas, previous$genotype, previous_snp, columns
    )
    same_matrix(previous_snp, snp)
    same_matrix(previous_peak$matrix, peak$matrix)
    stopifnot(identical(previous_peak$overlap, peak$overlap))
    same_matrix(previous_aligned, aligned)
    previous_snp_ld <- baseline$compute_snp_ld_covariance(previous_aligned, previous$gwas, columns)
    previous_peak_ld <- baseline$build_peak_ld_model(
      previous$gwas, previous_aligned, previous_peak$matrix, previous_peak$overlap
    )
    same_matrix(previous_snp_ld, snp_ld)
    same_numeric(baseline$compute_inverse_sqrt(previous_snp_ld), modified$compute_inverse_sqrt(snp_ld))
    for (field in c("response", "sigma", "covariance_q", "inverse_sqrt")) {
      same_numeric(previous_peak_ld[[field]], peak_ld[[field]])
    }
  }

  history_changed <- current$gwas
  history_changed$EPIC2_input_pos <- 100L
  history_changed$EPIC2_input_chr <- "2"
  history_changed$EPIC2_input_genome <- if (build == "hg19") "hg38" else "hg19"
  same_matrix(snp, modified$build_snp_accessibility(
    atac, history_changed, columns, workers = 1L, chunks = 2L, genome_build = build
  ))
  history_peak <- modified$build_peak_accessibility(atac, history_changed, columns, genome_build = build)
  same_matrix(peak$matrix, history_peak$matrix)
  stopifnot(identical(peak$overlap, history_peak$overlap))

  other_build <- if (build == "hg19") "hg38" else "hg19"
  conflicting_gwas <- current$gwas
  attr(conflicting_gwas, "genome") <- other_build
  expect_error_message(modified$build_peak_accessibility(
    atac, conflicting_gwas, columns, genome_build = build
  ), "genome metadata reports")
  conflicting_gwas <- current$gwas
  conflicting_gwas$EPIC2_target_genome <- other_build
  expect_error_message(modified$build_snp_accessibility(
    atac, conflicting_gwas, columns, workers = 1L, chunks = 2L, genome_build = build
  ), "EPIC2_target_genome")
  conflicting_atac <- atac
  conflicting_assay <- atac[["ATAC"]]
  GenomeInfoDb::genome(conflicting_assay@ranges) <- other_build
  conflicting_atac[["ATAC"]] <- conflicting_assay
  expect_error_message(modified$build_peak_accessibility(
    conflicting_atac, current$gwas, columns, genome_build = build
  ), "stored Seqinfo genome")
  message("Synthetic ", build, " accessibility and LD tests completed",
          if (!is.null(baseline)) " with historical numerical comparison." else ".")
}
run_synthetic_accessibility_test("hg19")
run_synthetic_accessibility_test("hg38")
