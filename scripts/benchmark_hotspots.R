#!/usr/bin/env Rscript

if (!requireNamespace("EPIC2", quietly = TRUE)) {
  stop("Install EPIC2 before running this benchmark.", call. = FALSE)
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Install Matrix before running this benchmark.", call. = FALSE)
}

set.seed(1000)
median_elapsed <- function(callback, repeats = 2L) {
  median(replicate(repeats, system.time(callback())[["elapsed"]]))
}

n_variants <- 3000L
n_samples <- 1000L
genotype <- matrix(
  sample(c(0, 1, 2, NA), n_variants * n_samples, replace = TRUE,
         prob = c(0.49, 0.38, 0.12, 0.01)),
  nrow = n_variants
)
rownames(genotype) <- paste0("rs", seq_len(n_variants))
colnames(genotype) <- paste0("sample", seq_len(n_samples))
gwas <- data.frame(
  rsid = rownames(genotype),
  chr = 1L,
  pos = seq.int(100000L, by = 1000L, length.out = n_variants),
  stringsAsFactors = FALSE
)
gwas$coord_id <- paste0("chr1:", gwas$pos)
rownames(genotype) <- gwas$coord_id

reference_scale_genotype <- function() {
  x <- t(as.matrix(genotype))
  for (column in seq_len(ncol(x))) {
    mean_value <- mean(x[, column], na.rm = TRUE)
    x[is.na(x[, column]), column] <- mean_value
    x[, column] <- (x[, column] - mean_value) / stats::sd(x[, column])
  }
  x
}

reference_scaled <- reference_scale_genotype()
optimized_scaled <- EPIC2:::scale_genotype_for_ld(genotype)
stopifnot(isTRUE(all.equal(reference_scaled, optimized_scaled, tolerance = 1e-12)))
scale_reference_seconds <- median_elapsed(reference_scale_genotype)
scale_optimized_seconds <- median_elapsed(function() {
  EPIC2:::scale_genotype_for_ld(genotype)
})

ld_scalar <- EPIC2::compute_snp_ld_covariance(
  genotype, gwas, window_bp = 500000L, correlation_block_size = 1L
)
ld_block <- EPIC2::compute_snp_ld_covariance(
  genotype, gwas, window_bp = 500000L, correlation_block_size = 256L
)
stopifnot(max(abs(ld_scalar - ld_block)) < 1e-12)
ld_scalar_seconds <- median_elapsed(function() {
  EPIC2::compute_snp_ld_covariance(
    genotype, gwas, window_bp = 500000L, correlation_block_size = 1L
  )
})
ld_block_seconds <- median_elapsed(function() {
  EPIC2::compute_snp_ld_covariance(
    genotype, gwas, window_bp = 500000L, correlation_block_size = 256L
  )
})

result <- data.frame(
  hotspot = c("genotype scaling", "local SNP LD"),
  reference_seconds = c(scale_reference_seconds, ld_scalar_seconds),
  optimized_seconds = c(scale_optimized_seconds, ld_block_seconds),
  speedup = c(
    scale_reference_seconds / scale_optimized_seconds,
    ld_scalar_seconds / ld_block_seconds
  )
)
print(result, row.names = FALSE, digits = 3)
