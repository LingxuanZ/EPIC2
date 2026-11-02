library(EPIC2)

expect_guard_error <- function(expression, text) {
  error <- tryCatch({ force(expression); NA_character_ }, error = conditionMessage)
  stopifnot(!is.na(error), all(vapply(text, grepl, logical(1), x = error, fixed = TRUE)))
  invisible(error)
}

# Small dimensions test the threshold without requiring a large allocation.
x <- matrix(c(1, 0.2, 0.1, 0.2, 1, 0.3, 0.1, 0.3, 1), 3L)
dimnames(x) <- list(paste0("s", 1:3), paste0("s", 1:3))
at_limit <- compute_inverse_sqrt(x, max_dense_block_size = 3L)
unlimited <- compute_inverse_sqrt(x, max_dense_block_size = Inf)
eig <- eigen((x + t(x)) / 2, symmetric = TRUE)
reference <- eig$vectors %*% diag(1 / sqrt(pmax(eig$values, 1e-8) + 0.01), 3L) %*%
  t(eig$vectors)
dimnames(reference) <- dimnames(x)
stopifnot(identical(at_limit, unlimited), identical(at_limit, reference))
expect_guard_error(compute_inverse_sqrt(x, max_dense_block_size = 2L),
                   c("3 x 3", "one double matrix", "planning estimate", "Inf"))

# A conversion sentinel proves rejection occurs before as.matrix(), not after.
as.matrix.epic2_guard_probe <- function(x, ...) stop("Unexpected dense conversion")
probe <- x
class(probe) <- c("epic2_guard_probe", class(probe))
expect_guard_error(compute_inverse_sqrt(probe, max_dense_block_size = 2L),
                   "exceeds max_dense_block_size")

for (invalid in list(0, -1, 1.5, NA_real_, NaN, -Inf, "5000", TRUE, c(2, 3), NULL)) {
  expect_guard_error(compute_inverse_sqrt(x, max_dense_block_size = invalid),
                     "positive integer or Inf")
}

# Configuration validation permits explicit unlimited jobs and rejects typos.
config <- EPIC2:::EPIC2_defaults()
config$genome$build <- "hg19"
config$atac$input <- "unused.rds"
config$gwas$traits <- list(TEST = list(summary = "unused.tsv", genotype = "unused.rds"))
stopifnot(config$ld$max_dense_block_size == 5000L)
validate_EPIC2_config(config)
config$ld$max_dense_block_size <- yaml::yaml.load("limit: .inf")$limit
stopifnot(identical(config$ld$max_dense_block_size, Inf))
validate_EPIC2_config(config)
for (invalid in list(0, -Inf, 1.5, NA_real_, "unlimited", TRUE, c(1, 2), NULL)) {
  config$ld$max_dense_block_size <- invalid
  expect_guard_error(validate_EPIC2_config(config), "ld.max_dense_block_size")
}

# Within-peak rejection precedes construction of its dense SNP correlation.
# No genome parsing or optional package work is reached in this failing case.
peak_function <- EPIC2:::build_peak_ld_model_impl
environment(peak_function) <- list2env(list(
  require_namespaces = function(...) invisible(TRUE),
  ld_from_scaled = function(...) stop("Unexpected within-peak dense allocation")
), parent = environment(peak_function))
genotype <- rbind(c(0, 1, 2, 1), c(2, 0, 1, 1), c(1, 2, 0, 1))
dimnames(genotype) <- list(rownames(x), paste0("sample", 1:4))
gwas <- data.frame(coord_id = rownames(x), Z = c(1, 2, 3))
peak_matrix <- matrix(1, 1L, 2L, dimnames = list("chr1-90-120", c("c1", "c2")))
overlap <- data.frame(snp = rownames(x), peak = rownames(peak_matrix))
expect_guard_error(
  peak_function(gwas, genotype, peak_matrix, overlap, max_dense_block_size = 2L),
  c("Within-peak SNP LD", "3 x 3", "exceeds max_dense_block_size")
)

if (requireNamespace("igraph", quietly = TRUE)) {
  # A chain of local nonzero links is one component despite sparse storage.
  chain <- Matrix::sparseMatrix(i = c(1:3, 1, 2, 2, 3),
                                j = c(1:3, 2, 1, 3, 2),
                                x = c(rep(1, 3), rep(0.2, 4)), dims = c(3L, 3L))
  small <- matrix(c(1, 0.1, 0.1, 1), 2L)
  sigma <- EPIC2:::as_general_csparse(Matrix::bdiag(small, chain))
  dimnames(sigma) <- list(paste0("v", 1:5), paste0("v", 1:5))
  component_function <- EPIC2:::inverse_sqrt_components
  environment(component_function) <- list2env(list(
    split_sparse_components = function(...) list(small = paste0("v", 1:2),
                                                 large = paste0("v", 3:5)),
    compute_inverse_sqrt = function(...) stop("Unexpected first eigen block")
  ), parent = environment(component_function))
  expect_guard_error(component_function(sigma, 0.01, 1e-8, max_dense_block_size = 2L),
                     c("Largest LD connected component", "3 x 3"))
  # The total matrix dimension may exceed the limit when each block fits.
  accepted <- EPIC2:::inverse_sqrt_components(sigma, 0.01, 1e-8, max_dense_block_size = 3L)
  original <- EPIC2:::inverse_sqrt_components(sigma, 0.01, 1e-8, max_dense_block_size = Inf)
  stopifnot(identical(accepted, original))

  peak_dependencies <- c("Signac", "GenomicRanges")
  if (all(vapply(peak_dependencies, requireNamespace, logical(1), quietly = TRUE))) {
    # The cache limit changes storage strategy, with the same peak covariance.
    peak_ids <- c("chr1-90-120", "chr1-190-220")
    two_peaks <- matrix(1, 2, 2, dimnames = list(peak_ids, c("c1", "c2")))
    two_overlap <- data.frame(snp = rownames(genotype), peak = peak_ids[c(1, 1, 2)])
    cached <- build_peak_ld_model(gwas, genotype, two_peaks, two_overlap,
                                  max_dense_block_size = Inf)
    uncached <- build_peak_ld_model(gwas, genotype, two_peaks, two_overlap,
                                    max_dense_block_size = 2L)
    stopifnot(isTRUE(all.equal(cached, uncached, tolerance = 1e-12)))
    # A limit that allows each within-peak inverse still rejects a larger
    # connected peak-response covariance before decomposition.
    singleton_overlap <- data.frame(snp = rownames(genotype)[1:2], peak = peak_ids)
    expect_guard_error(
      build_peak_ld_model(gwas, genotype, two_peaks, singleton_overlap,
                          max_dense_block_size = 1L),
      c("Largest LD connected component", "2 x 2")
    )
  } else message("Peak-cache equivalence test requires Signac and GenomicRanges.")
} else message("Connected-component guard test requires igraph.")

message("LD dense-memory guard tests passed.")
