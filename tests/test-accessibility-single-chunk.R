# Regression tests use actual indexed fragments, so successful execution checks
# the complete FeatureMatrix path as well as the one-chunk split.
source_directory <- Sys.getenv("EPIC2_MODIFIED_R_DIR", "R")
if (dir.exists(source_directory)) {
  epic2 <- new.env(parent = globalenv())
  for (path in sort(list.files(source_directory, pattern = "\\.R$", full.names = TRUE))) {
    sys.source(path, envir = epic2)
  }
} else {
  epic2 <- asNamespace("EPIC2")
}

required <- c("Signac", "Seurat", "Rsamtools", "Matrix", "GenomicRanges", "IRanges",
              "GenomeInfoDb", "BiocParallel")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Single-chunk accessibility test requires: ", paste(missing, collapse = ", "))
}

run_single_chunk_regression <- function() {
  directory <- tempfile("EPIC2-single-chunk-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  cells <- paste0("cell", seq_len(3L))
  peaks <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(start = c(999900L, 1000900L), end = c(1000100L, 1001100L))
  )
  GenomeInfoDb::genome(peaks) <- "hg19"
  fragments_table <- data.frame(
    chr = "chr1", start = c(999995L, 999995L, 1000995L, 1000995L),
    end = c(1000005L, 1000005L, 1001005L, 1001005L),
    cell = c("cell1", "cell2", "cell2", "cell3"), count = 1L
  )
  fragment_path <- file.path(directory, "fragments.tsv")
  utils::write.table(fragments_table, fragment_path, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  fragment_gz <- Rsamtools::bgzip(fragment_path, overwrite = TRUE)
  Rsamtools::indexTabix(fragment_gz, format = "bed")
  fragments <- Signac::CreateFragmentObject(
    fragment_gz, cells = setNames(cells, cells),
    validate.fragments = TRUE, verbose = FALSE
  )
  counts <- Signac::FeatureMatrix(fragments, features = peaks, cells = cells, verbose = FALSE)
  chromatin <- Signac::CreateChromatinAssay(
    counts, ranges = peaks, fragments = fragments, min.cells = 0, min.features = 0,
    validate.fragments = FALSE
  )
  atac <- Seurat::CreateSeuratObject(counts = chromatin, assay = "ATAC")
  gwas <- data.frame(
    rsid = c("rs2", "rs1"), chr = 1L, pos = c(1001000L, 1000000L),
    coord_id = c("chr1-1001000-1001000", "chr1-1000000-1000000")
  )
  expected <- rbind(c(0, 1, 1), c(1, 1, 0))
  dimnames(expected) <- list(gwas$coord_id, cells)
  assert_expected <- function(observed, expected) {
    stopifnot(
      identical(dim(observed), dim(expected)),
      identical(as.vector(rownames(observed)), as.vector(rownames(expected))),
      identical(as.vector(colnames(observed)), as.vector(colnames(expected))),
      isTRUE(all.equal(as.matrix(observed), expected, check.attributes = FALSE)),
      identical(attr(observed, "genome_build"), "hg19")
    )
  }
  # One retained SNP must work with the default chunks=20.
  single <- epic2$build_snp_accessibility(atac, gwas[1L, , drop = FALSE], genome_build = "hg19")
  assert_expected(single, expected[1L, , drop = FALSE])
  # Explicit chunks=1 must work with multiple retained SNPs and preserve GWAS order.
  one_chunk <- epic2$build_snp_accessibility(atac, gwas, chunks = 1L, genome_build = "hg19")
  assert_expected(one_chunk, expected)
  two_chunks <- epic2$build_snp_accessibility(atac, gwas, chunks = 2L, genome_build = "hg19")
  assert_expected(two_chunks, expected)
}

run_single_chunk_regression()
message("Single-SNP and single-chunk accessibility regression tests passed.")
