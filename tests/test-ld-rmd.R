library(EPIC2)

# Compare named mathematical results, not incidental sparse storage or block order.
expect_named_matrix_equal <- function(actual, expected, tolerance = 1e-8) {
  stopifnot(
    setequal(rownames(actual), rownames(expected)),
    setequal(colnames(actual), colnames(expected))
  )
  actual <- as.matrix(actual[rownames(expected), colnames(expected), drop = FALSE])
  stopifnot(max(abs(actual - as.matrix(expected))) < tolerance)
}

# Notebook reference on complete genotypes. Missing-genotype policies are not
# equivalent: the package mean-imputes, whereas the notebook retains NA values.
rmd_scale_complete_genotype <- function(genotype) {
  stopifnot(all(is.finite(genotype)))
  x <- t(as.matrix(genotype))
  standard_deviation <- apply(x, 2, stats::sd)
  keep <- is.finite(standard_deviation) & standard_deviation > 0
  scaled <- as.matrix(scale(x[, keep, drop = FALSE]))
  colnames(scaled) <- rownames(genotype)[keep]
  scaled
}

rmd_inverse_sqrt <- function(x, ridge = 0.01, eigen_tolerance = 1e-8) {
  x <- as.matrix(x)
  decomposition <- eigen((x + t(x)) / 2, symmetric = TRUE)
  values <- pmax(decomposition$values, eigen_tolerance) + ridge
  result <- decomposition$vectors %*%
    diag(1 / sqrt(values), nrow = length(values)) %*%
    t(decomposition$vectors)
  dimnames(result) <- dimnames(x)
  result
}

# Independent scalar-focal reference from the active 5 Mb SNP notebook cell.
# A smaller window makes the inclusive boundary easy to exercise in a fixture.
rmd_snp_covariance <- function(genotype, gwas, window_bp, factor_two = FALSE) {
  genotype <- genotype[match(gwas$coord_id, rownames(genotype)), , drop = FALSE]
  ord <- order(as.integer(gwas$chr), gwas$pos)
  gwas <- gwas[ord, , drop = FALSE]
  genotype <- genotype[ord, , drop = FALSE]
  scaled <- rmd_scale_complete_genotype(genotype)
  gwas <- gwas[gwas$coord_id %in% colnames(scaled), , drop = FALSE]
  scaled <- scaled[, gwas$coord_id, drop = FALSE]
  sigma <- matrix(
    0, nrow(gwas), nrow(gwas),
    dimnames = list(gwas$coord_id, gwas$coord_id)
  )
  for (index in split(seq_len(nrow(gwas)), gwas$chr)) {
    positions <- gwas$pos[index]
    chromosome <- scaled[, gwas$coord_id[index], drop = FALSE]
    for (a in seq_along(index)) {
      right <- findInterval(positions[a] + window_bp, positions)
      js <- seq.int(a, right)
      correlation <- as.numeric(crossprod(
        chromosome[, a, drop = FALSE],
        chromosome[, js, drop = FALSE]
      ) / (nrow(chromosome) - 1))
      value <- correlation^2 * if (factor_two) 2 else 1
      sigma[index[a], index[js]] <- value
      sigma[index[js], index[a]] <- value
    }
  }
  diag(sigma) <- if (factor_two) 2 else 1
  sigma
}

ld_ids <- c(
  "chr1-100-101", "chr1-600-601", "chr1-601-602",
  "chr1-2000-2001", "chr2-100-101", "chr1-250-251"
)
ld_genotype <- rbind(
  c(0, 0, 1, 1, 2, 2, 1, 0),
  c(0, 0, 1, 1, 2, 2, 1, 0),
  c(2, 2, 1, 1, 0, 0, 1, 2),
  c(0, 1, 0, 2, 1, 0, 2, 1),
  c(0, 0, 1, 1, 2, 2, 1, 0),
  rep(1, 8)
)
dimnames(ld_genotype) <- list(ld_ids, paste0("sample", seq_len(8)))
ld_gwas <- data.frame(
  coord_id = ld_ids,
  chr = c(1, 1, 1, 1, 2, 1),
  pos = c(100, 600, 601, 2000, 100, 250)
)
ld_genotype <- ld_genotype[c(6, 3, 1, 5, 2, 4), , drop = FALSE]
ld_gwas <- ld_gwas[c(5, 2, 6, 4, 1, 3), , drop = FALSE]

expected_snp <- rmd_snp_covariance(ld_genotype, ld_gwas, window_bp = 500)
for (block_size in c(1L, 2L, 256L)) {
  actual_snp <- compute_snp_ld_covariance(
    ld_genotype, ld_gwas,
    window_bp = 500, correlation_block_size = block_size
  )
  expect_named_matrix_equal(actual_snp, expected_snp, tolerance = 1e-12)
  stopifnot(
    !ld_ids[6] %in% rownames(actual_snp),
    abs(actual_snp[ld_ids[1], ld_ids[2]] - 1) < 1e-12,
    actual_snp[ld_ids[1], ld_ids[3]] == 0,
    actual_snp[ld_ids[1], ld_ids[5]] == 0
  )
  actual_inverse <- EPIC2:::inverse_sqrt_components(actual_snp, 0.01, 1e-8)
  expect_named_matrix_equal(actual_inverse$matrix, rmd_inverse_sqrt(expected_snp))
}
actual_factor_two <- compute_snp_ld_covariance(
  ld_genotype, ld_gwas, window_bp = 500, factor_two = TRUE
)
expect_named_matrix_equal(actual_factor_two, 2 * expected_snp, tolerance = 1e-12)

# Explicit stored zeros must not become graph edges. The notebook's sparse
# logical conversion can retain FALSE entries; reproducing that bug is not needed
# for equality of the named covariance/inverse-square-root results.
isolated <- Matrix::sparseMatrix(
  i = c(1, 1, 2, 2, 3), j = c(1, 2, 1, 2, 3),
  x = c(1, 0, 0, 1, 1),
  dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
)
isolated_inverse <- EPIC2:::inverse_sqrt_components(isolated, 0.01, 1e-8)
stopifnot(all(lengths(isolated_inverse$indices) == 1L))
expect_named_matrix_equal(isolated_inverse$matrix, rmd_inverse_sqrt(isolated))
message("LD notebook regression: SNP covariance and inverse-square-root checks passed.")

# Stage-level checkpoint ordering is testable without Bioconductor peak parsing.
# Clone the stage function and fail exactly at entry to peak LD, after SNP saves.
local({
  stage_dir <- tempfile("EPIC2-ld-stage-checkpoint-")
  dir.create(stage_dir)
  on.exit(unlink(stage_dir, recursive = TRUE), add = TRUE)
  dir.create(file.path(stage_dir, "02_gwas"))
  dir.create(file.path(stage_dir, "02_accessibility"))

  rsid_lookup <- setNames(paste0("rs", seq_along(ld_ids)), ld_ids)
  stage_gwas <- ld_gwas
  stage_gwas$rsid <- unname(rsid_lookup[stage_gwas$coord_id])
  stage_gwas$Z <- 1
  stage_genotype <- ld_genotype
  rownames(stage_genotype) <- unname(rsid_lookup[rownames(stage_genotype)])
  stage_snps <- Matrix::Matrix(
    matrix(1, length(ld_ids), 2, dimnames = list(ld_ids, c("cell1", "cell2"))),
    sparse = TRUE
  )
  stage_peaks <- Matrix::Matrix(
    matrix(1, 1, 2, dimnames = list("chr1-90-120", c("cell1", "cell2"))),
    sparse = TRUE
  )
  saveRDS(stage_gwas, file.path(stage_dir, "02_gwas", "HDL.gwas.rds"))
  saveRDS(stage_genotype, file.path(stage_dir, "02_gwas", "HDL.genotype.rds"))
  saveRDS(stage_snps, file.path(stage_dir, "02_accessibility", "HDL.snp_by_cell.rds"))
  saveRDS(stage_peaks, file.path(stage_dir, "02_accessibility", "HDL.peak_by_cell.rds"))
  utils::write.table(
    data.frame(snp = ld_ids[1], peak = "chr1-90-120"),
    file.path(stage_dir, "02_accessibility", "HDL.snp_peak_overlap.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  config <- EPIC2:::EPIC2_defaults()
  config$output_dir <- stage_dir
  config$project$output_dir <- stage_dir
  config$gwas$traits <- list(HDL = list())
  config$ld$window_bp <- 500

  stage_function <- run_ld_stage
  environment(stage_function) <- list2env(
    list(build_peak_ld_model_impl = function(...) stop("Simulated peak-stage failure")),
    parent = environment(run_ld_stage)
  )
  failure <- tryCatch({
    stage_function(config, traits = "HDL")
    NA_character_
  }, error = conditionMessage)
  prefix <- file.path(stage_dir, "03_ld", "HDL")
  stopifnot(
    identical(failure, "Simulated peak-stage failure"),
    file.exists(paste0(prefix, ".snp_sigma.rds")),
    file.exists(paste0(prefix, ".snp_inverse_sqrt.rds")),
    file.exists(paste0(prefix, ".snp_blocks.rds")),
    file.exists(paste0(prefix, ".ld_runtime.tsv")),
    !file.exists(paste0(prefix, ".peak_ld_model.rds"))
  )
  expect_named_matrix_equal(readRDS(paste0(prefix, ".snp_sigma.rds")), expected_snp)
  expect_named_matrix_equal(
    readRDS(paste0(prefix, ".snp_inverse_sqrt.rds")), rmd_inverse_sqrt(expected_snp)
  )
  runtime <- utils::read.delim(paste0(prefix, ".ld_runtime.tsv"), stringsAsFactors = FALSE)
  stopifnot(identical(runtime$value[runtime$setting == "status"], "RUNNING"))
})
message("LD notebook regression: SNP files persist after simulated peak-stage failure.")

peak_dependencies <- c("Signac", "GenomicRanges")
if (!all(vapply(peak_dependencies, requireNamespace, logical(1), quietly = TRUE))) {
  message("Skipping peak LD notebook regression: Signac/GenomicRanges are unavailable.")
} else {
  # Independent dense reference for the finite-input notebook peak equations.
  # Peak pairs use start(next) <= end(current) + window, not start-to-start distance.
  rmd_peak_reference <- function(genotype, gwas, peaks, overlap,
                                  window_bp = 500, ridge = 0.01,
                                  eigen_tolerance = 1e-8) {
    scaled <- rmd_scale_complete_genotype(genotype)
    z_lookup <- setNames(gwas$Z, gwas$coord_id)
    overlap <- unique(overlap[
      overlap$snp %in% colnames(scaled) & overlap$peak %in% peaks,
      c("snp", "peak"), drop = FALSE
    ])
    snp_sets <- lapply(split(overlap$snp, overlap$peak), unique)
    peak_ids <- names(snp_sets)
    within <- lapply(snp_sets, function(ids) {
      r <- crossprod(scaled[, ids, drop = FALSE]) / (nrow(scaled) - 1)
      if (length(ids) == 1L) {
        matrix(1 / sqrt(r[1, 1] + ridge), 1, 1, dimnames = list(ids, ids))
      } else {
        rmd_inverse_sqrt(r, ridge, eigen_tolerance)
      }
    })
    k <- lengths(snp_sets)
    response <- vapply(peak_ids, function(peak) {
      sum(as.numeric(within[[peak]] %*% z_lookup[snp_sets[[peak]]])^2) / k[[peak]]
    }, numeric(1))
    pieces <- strsplit(peak_ids, "-", fixed = TRUE)
    chromosome <- vapply(pieces, `[[`, character(1), 1L)
    start <- as.integer(vapply(pieces, `[[`, character(1), 2L))
    end <- as.integer(vapply(pieces, `[[`, character(1), 3L))
    covariance_q <- matrix(
      0, length(peak_ids), length(peak_ids), dimnames = list(peak_ids, peak_ids)
    )
    for (index in split(order(start, end), chromosome[order(start, end)])) {
      for (a in seq_along(index)) {
        i <- index[a]
        for (b in seq.int(a, length(index))) {
          j <- index[b]
          if (start[j] > end[i] + window_bp) next
          snps1 <- snp_sets[[i]]
          snps2 <- snp_sets[[j]]
          r12 <- crossprod(
            scaled[, snps1, drop = FALSE], scaled[, snps2, drop = FALSE]
          ) / (nrow(scaled) - 1)
          transformed <- within[[i]] %*% r12 %*% within[[j]]
          covariance_q[i, j] <- 2 * sum(transformed * transformed)
          covariance_q[j, i] <- covariance_q[i, j]
        }
      }
    }
    sigma <- covariance_q / outer(k, k)
    list(
      response = response, covariance_q = covariance_q, sigma = sigma,
      inverse_sqrt = rmd_inverse_sqrt(sigma, ridge, eigen_tolerance)
    )
  }

  peak_variant_ids <- c(
    "chr1-100-101", "chr1-110-111", "chr1-600-601",
    "chr1-601-602", "chr2-100-101", "chr1-105-106"
  )
  peak_genotype <- rbind(
    c(0, 0, 1, 1, 2, 2, 1, 0),
    c(0, 0, 1, 1, 2, 2, 1, 0),
    c(0, 0, 1, 1, 2, 2, 1, 0),
    c(2, 2, 1, 1, 0, 0, 1, 2),
    c(0, 1, 0, 2, 1, 0, 2, 1),
    rep(1, 8)
  )
  dimnames(peak_genotype) <- list(peak_variant_ids, paste0("sample", seq_len(8)))
  peak_gwas <- data.frame(coord_id = peak_variant_ids, Z = c(2, 3, 1, -2, 0.5, 4))
  peak_ids <- c(
    "chr1-99-100", "chr1-90-120", "chr1-600-600",
    "chr1-601-601", "chr2-90-120"
  )
  peak_matrix <- Matrix::Matrix(
    matrix(1, 5, 3, dimnames = list(peak_ids[c(5, 4, 1, 2, 3)], paste0("cell", 1:3))),
    sparse = TRUE
  )
  peak_overlap <- data.frame(
    snp = peak_variant_ids[c(1, 1, 2, 3, 4, 5, 6, 1)],
    peak = peak_ids[c(1, 2, 2, 3, 4, 5, 2, 1)]
  )
  peak_genotype <- peak_genotype[c(6, 3, 1, 5, 2, 4), , drop = FALSE]
  peak_gwas <- peak_gwas[c(5, 2, 6, 4, 1, 3), , drop = FALSE]
  expected_peak <- rmd_peak_reference(
    peak_genotype, peak_gwas, rownames(peak_matrix), peak_overlap
  )

  for (cache_limit in c(1L, 15000L)) {
    actual_peak <- build_peak_ld_model(
      peak_gwas, peak_genotype, peak_matrix, peak_overlap,
      window_bp = 500, max_snps_for_chr_ld = cache_limit
    )
    stopifnot(
      setequal(names(actual_peak$response), names(expected_peak$response)),
      max(abs(actual_peak$response[names(expected_peak$response)] - expected_peak$response)) < 1e-8,
      actual_peak$peak_info[[peak_ids[2]]]$k == 2L,
      actual_peak$sigma[peak_ids[1], peak_ids[3]] > 0,
      actual_peak$sigma[peak_ids[1], peak_ids[4]] == 0,
      actual_peak$sigma[peak_ids[1], peak_ids[5]] == 0
    )
    for (component in c("covariance_q", "sigma", "inverse_sqrt")) {
      expect_named_matrix_equal(actual_peak[[component]], expected_peak[[component]])
    }
  }

  # Forty-eight mutually overlapping peaks generate 1,176 accepted upper-triangle
  # pairs, exceeding the initial 1,024-entry capacity and exercising array growth.
  growth_ids <- paste0("chr1-90-", seq.int(100, 147))
  growth_matrix <- Matrix::Matrix(
    matrix(1, 48, 2, dimnames = list(growth_ids, c("cell1", "cell2"))),
    sparse = TRUE
  )
  growth_overlap <- data.frame(snp = peak_variant_ids[1], peak = growth_ids)
  growth_model <- build_peak_ld_model(
    peak_gwas, peak_genotype, growth_matrix, growth_overlap, window_bp = 500
  )
  expected_growth_sigma <- matrix(
    2 / (1 + 0.01)^2, 48, 48, dimnames = list(growth_ids, growth_ids)
  )
  stopifnot(
    48L * 49L / 2L > 1024L,
    Matrix::nnzero(growth_model$sigma) == 48L * 48L,
    max(abs(growth_model$response - 2^2 / (1 + 0.01))) < 1e-12
  )
  expect_named_matrix_equal(growth_model$sigma, expected_growth_sigma, tolerance = 1e-12)
  expect_named_matrix_equal(
    growth_model$inverse_sqrt, rmd_inverse_sqrt(expected_growth_sigma)
  )

  # A simulated failure in the final inverse step must leave earlier checkpoints.
  # A cloned closure overrides only this call; package namespace state is untouched.
  local({
    checkpoint_dir <- tempfile("EPIC2-ld-checkpoint-")
    dir.create(checkpoint_dir)
    on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
    original <- EPIC2:::build_peak_ld_model_impl
    peak_function <- original
    environment(peak_function) <- list2env(
      list(inverse_sqrt_components = function(...) stop("Simulated peak inverse failure")),
      parent = environment(original)
    )
    prefix <- file.path(checkpoint_dir, "HDL")
    failure <- tryCatch({
      peak_function(
        peak_gwas, peak_genotype, peak_matrix, peak_overlap,
        window_bp = 500, checkpoint_prefix = prefix
      )
      NA_character_
    }, error = conditionMessage)
    stopifnot(
      identical(failure, "Simulated peak inverse failure"),
      file.exists(paste0(prefix, ".peak_response.rds")),
      file.exists(paste0(prefix, ".peak_covariance.rds"))
    )
  })
  message("LD notebook regression: peak equations, cache/fallback and early checkpoints passed.")
}
