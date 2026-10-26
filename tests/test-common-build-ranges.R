# Run from the source root, or against the installed package under R CMD check.
validation_env <- new.env(parent = globalenv())
if (file.exists(file.path("R", "genome.R"))) {
  sys.source(file.path("R", "utils.R"), envir = validation_env)
  sys.source(file.path("R", "genome.R"), envir = validation_env)
} else {
  namespace <- asNamespace("EPIC2")
  validation_env <- new.env(parent = namespace)
  for (name in ls(namespace, all.names = TRUE)) {
    value <- get(name, envir = namespace, inherits = FALSE)
    if (is.function(value) && identical(environment(value), namespace)) {
      environment(value) <- validation_env
      assign(name, value, envir = validation_env)
    }
  }
}
expect_error_message <- function(expression, pattern) {
  message <- tryCatch({ force(expression); NA_character_ }, error = conditionMessage)
  stopifnot(!is.na(message), grepl(pattern, message, ignore.case = TRUE))
}

common_build <- validation_env$assert_common_genome_build
stopifnot(
  identical(common_build("GRCh37", list("atac.genome" = "hg19", "gwas.genome" = NULL)), "hg19"),
  identical(common_build("hg38", list("gwas.traits.HDL.genome" = "GRCh38")), "hg38")
)
expect_error_message(common_build(), "genome.build is required")
expect_error_message(common_build(NULL), "genome.build is required")
expect_error_message(common_build(""), "invalid genome.build")
expect_error_message(common_build("GRCh19"), "Unsupported genome build")
expect_error_message(common_build("hg19", list("atac.genome" = "hg38")), "atac.genome.*hg38.*hg19")
expect_error_message(common_build("hg19", list("gwas.traits.HDL.genome" = NA_character_)),
                     "invalid gwas.traits.HDL.genome")

required <- c("GenomicRanges", "IRanges", "S4Vectors", "GenomeInfoDb")
stopifnot(all(vapply(required, requireNamespace, logical(1), quietly = TRUE)))
ranges <- GenomicRanges::GRanges(
  seqnames = c("2", "1", "chr1", "chr2"),
  ranges = IRanges::IRanges(start = c(300L, 100L, 100L, 50L), end = c(310L, 110L, 110L, 60L))
)
names(ranges) <- c("last_chromosome", "duplicate_a", "duplicate_b", "earlier_position")
S4Vectors::mcols(ranges)$input_row <- seq_along(ranges)
old_ranges <- validation_env$lift_genomic_ranges(ranges, "hg19", "hg19")
old_ranges_hg38 <- validation_env$lift_genomic_ranges(ranges, "hg38", "hg38")

gwas <- data.frame(
  rsid = paste0("rs", 1:4), chr = c("2", "1", "chr1", "chr2"),
  pos = c(300L, 100L, 100L, 50L), A1 = c("A", "C", "G", "T"),
  Beta = c(0.3, 0.1, -0.2, 0.5), build = "GRCh37",
  EPIC2_input_pos = 999L, stringsAsFactors = FALSE
)
columns <- list(chr = "chr", pos = "pos")
old_gwas <- validation_env$harmonize_gwas_coordinates(
  gwas, columns, "hg19", "hg19", genome_column = "build"
)

# Any accidental delegation to coordinate conversion must make these tests fail.
validation_env$lift_genomic_ranges <- function(...) stop("Unexpected liftOver call")
validation_env$harmonize_gwas_coordinates <- function(...) stop("Unexpected harmonization call")
new_ranges <- validation_env$validate_genomic_ranges(ranges, "GRCh37")
new_ranges_hg38 <- validation_env$validate_genomic_ranges(ranges, "GRCh38")
stopifnot(
  identical(new_ranges$ranges, old_ranges$ranges),
  identical(new_ranges$input_index, c(1L, 4L)),
  identical(new_ranges$input_index, old_ranges$input_index),
  identical(new_ranges_hg38$ranges, old_ranges_hg38$ranges),
  new_ranges$qc$n_duplicate_target == 2L,
  identical(new_ranges$qc$operation, "validation_only"),
  identical(new_ranges$qc$conversion_performed, FALSE),
  is.na(new_ranges$qc$chain_file), is.na(new_ranges$qc$chain_md5)
)
expect_error_message(validation_env$validate_genomic_ranges(ranges, "hg19", "error"), "duplicated")
expect_error_message(validation_env$validate_genomic_ranges(ranges[2:3], "hg19"), "no unique ranges")
bad_names <- ranges
names(bad_names)[2L] <- names(bad_names)[1L]
expect_error_message(validation_env$validate_genomic_ranges(bad_names, "hg19"), "unique.*range names")
out_of_bounds <- GenomicRanges::GRanges("chr1", IRanges::IRanges(249000000L, 249000000L))
expect_error_message(validation_env$validate_genomic_ranges(out_of_bounds, "hg38"), "out-of-bounds")
stored_build <- ranges
GenomeInfoDb::genome(stored_build) <- "hg38"
expect_error_message(validation_env$validate_genomic_ranges(stored_build, "hg19"), "stored Seqinfo genome")

new_gwas <- validation_env$validate_gwas_coordinates(gwas, columns, "GRCh37", genome_column = "build")
comparable_gwas <- new_gwas$gwas
attr(comparable_gwas, "genome") <- NULL
stopifnot(
  identical(comparable_gwas, old_gwas$gwas),
  identical(new_gwas$input_index, old_gwas$input_index),
  identical(new_gwas$gwas$rsid, c("rs1", "rs4")),
  identical(new_gwas$gwas$EPIC2_input_pos, c(300L, 50L)),
  identical(new_gwas$gwas$A1, gwas$A1[new_gwas$input_index]),
  identical(attr(new_gwas$gwas, "genome"), "hg19"),
  identical(new_gwas$genome_qc$genome_build, "hg19"),
  identical(new_gwas$genome_qc$conversion_performed, FALSE),
  new_gwas$genome_qc$n_duplicate_target == 2L,
  identical(new_gwas$metadata_qc$metadata_status, "matched")
)
expect_error_message(
  validation_env$validate_gwas_coordinates(gwas, columns, "hg19", duplicate_action = "error"),
  "duplicated"
)
expect_error_message(
  validation_env$validate_gwas_coordinates(gwas, columns, "hg38", genome_column = "build"),
  "genome metadata reports hg19"
)
gwas_hg38 <- gwas
gwas_hg38$build <- "GRCh38"
stopifnot(identical(
  attr(validation_env$validate_gwas_coordinates(gwas_hg38, columns, "hg38", genome_column = "build")$gwas,
       "genome"), "hg38"
))
message("Common-build range and GWAS validation tests completed.")
