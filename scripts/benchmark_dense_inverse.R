#!/usr/bin/env Rscript

# Run after installing EPIC2:
# Rscript scripts/benchmark_dense_inverse.R [output.tsv] [250,500,1000]
# Memory columns are planning estimates, not measured peak RSS.
if (!requireNamespace("EPIC2", quietly = TRUE) ||
    !requireNamespace("Matrix", quietly = TRUE)) {
  stop("Install EPIC2 and Matrix before running this benchmark.", call. = FALSE)
}
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1L]] else "dense_inverse_benchmark.tsv"
sizes <- if (length(args) >= 2L) {
  suppressWarnings(as.numeric(strsplit(args[[2L]], ",", fixed = TRUE)[[1L]]))
} else c(250, 500, 1000)
if (!length(sizes) || any(!is.finite(sizes)) || any(sizes < 2) ||
    any(sizes != floor(sizes))) {
  stop("Sizes must be comma-separated integers of at least two.", call. = FALSE)
}
rows <- lapply(sizes, function(n) {
  # One connected, well-conditioned covariance component. Runtime is not a
  # claim about real reference-panel LD or a complete EPIC2 analysis.
  covariance <- Matrix::bandSparse(
    n, k = c(-1L, 0L, 1L),
    diagonals = list(rep(0.2, n - 1L), rep(1, n), rep(0.2, n - 1L))
  )
  ids <- paste0("feature", seq_len(n))
  dimnames(covariance) <- list(ids, ids)
  gc()
  elapsed <- system.time({
    inverse <- EPIC2::compute_inverse_sqrt(covariance, ridge = 0.01)
  })[["elapsed"]]
  probe <- sin(seq_len(n))
  error <- max(abs(as.numeric(inverse %*%
    ((covariance + Matrix::Diagonal(n, x = 0.01)) %*% (inverse %*% probe))) - probe))
  stopifnot(is.finite(error), error < 1e-8)
  data.frame(
    features = n, elapsed_seconds = elapsed,
    one_dense_matrix_gib = 8 * n^2 / 1024^3,
    six_dense_matrices_planning_gib = 6 * 8 * n^2 / 1024^3,
    probe_error = error,
    EPIC2_version = as.character(utils::packageVersion("EPIC2")),
    R_version = R.version.string,
    platform = R.version$platform,
    BLAS = unname(extSoftVersion()["BLAS"]),
    stringsAsFactors = FALSE
  )
})
result <- do.call(rbind, rows)
utils::write.table(result, output, sep = "\t", row.names = FALSE, quote = FALSE)
print(result[, c("features", "elapsed_seconds", "one_dense_matrix_gib", "probe_error")],
      row.names = FALSE)
cat("Saved ", output, ". Memory estimates exclude other live objects and are not peak-RSS measurements.\n", sep = "")
