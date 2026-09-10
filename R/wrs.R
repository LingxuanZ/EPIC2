#' Convert two-sided Wilcoxon results to enrichment-oriented one-sided scores
#'
#' @param result Data frame returned by `presto::wilcoxauc()`.
#' @param truncation_quantiles Lower and upper empirical winsorization quantiles.
#' @return Input data frame with one-sided P-value and Z-score columns.
#' @export
make_one_sided_wrs <- function(result, truncation_quantiles = c(0.02, 0.98)) {
  required <- c("pval", "auc")
  missing <- setdiff(required, colnames(result))
  if (length(missing)) stop("Missing Wilcoxon columns: ", paste(missing, collapse = ", "), call. = FALSE)
  result <- as.data.frame(result)
  p_two <- as.numeric(result$pval)
  auc <- as.numeric(result$auc)
  p_one <- ifelse(auc >= 0.5, p_two / 2, 1 - p_two / 2)
  p_for_z <- pmin(pmax(p_one, 1e-300), 1 - 1e-16)
  z <- stats::qnorm(p_for_z, lower.tail = FALSE)
  finite <- z[is.finite(z)]
  if (length(finite) == 0L) stop("No finite Wilcoxon Z-scores were generated.", call. = FALSE)
  limits <- stats::quantile(finite, truncation_quantiles, na.rm = TRUE, names = FALSE)
  result$pval_one_sided_enriched <- p_one
  result$z_one_sided_enriched <- z
  result$z_one_sided_enriched_trunc <- pmin(pmax(z, limits[[1L]]), limits[[2L]])
  result$z_trunc_lower <- limits[[1L]]
  result$z_trunc_upper <- limits[[2L]]
  result
}

run_wrs_matrix <- function(matrix, groups, truncation_quantiles) {
  require_namespaces("presto", "Wilcoxon rank-sum analysis")
  if (is.null(names(groups))) names(groups) <- colnames(matrix)
  groups <- groups[colnames(matrix)]
  keep <- !is.na(groups)
  if (length(unique(groups[keep])) < 2L) stop("Wilcoxon analysis requires at least two groups.", call. = FALSE)
  result <- presto::wilcoxauc(X = matrix[, keep, drop = FALSE], y = as.character(groups[keep]))
  make_one_sided_wrs(result, truncation_quantiles)
}

run_balanced_wrs_group <- function(matrix, groups, target, cells_per_side,
                                   repeats, seed, truncation_quantiles) {
  require_namespaces("presto", "balanced Wilcoxon analysis")
  if (is.null(names(groups))) names(groups) <- colnames(matrix)
  groups <- as.character(groups[colnames(matrix)])
  names(groups) <- colnames(matrix)
  inside <- names(groups)[!is.na(groups) & groups == target]
  outside <- names(groups)[!is.na(groups) & groups != target]
  if (length(inside) < cells_per_side || length(outside) < cells_per_side) return(NULL)
  set.seed(seed)
  results <- vector("list", repeats)
  for (repeat_index in seq_len(repeats)) {
    selected_inside <- sample(inside, cells_per_side, replace = FALSE)
    selected_outside <- sample(outside, cells_per_side, replace = FALSE)
    selected <- c(selected_inside, selected_outside)
    label <- c(rep("target", cells_per_side), rep("other", cells_per_side))
    result <- presto::wilcoxauc(X = matrix[, selected, drop = FALSE], y = label)
    result <- result[result$group == "target", , drop = FALSE]
    result <- make_one_sided_wrs(result, truncation_quantiles)
    result$target_group <- target
    result[["repeat"]] <- repeat_index
    result$cells_per_side <- cells_per_side
    results[[repeat_index]] <- result
  }
  do.call(rbind, results)
}

summarize_balanced_wrs <- function(result, auc_cutoff = 0.55) {
  if (is.null(result) || nrow(result) == 0L) return(data.frame())
  split_result <- split(result, interaction(result$target_group, result$feature, drop = TRUE))
  rows <- lapply(split_result, function(x) data.frame(
    group = x$target_group[[1L]],
    feature = x$feature[[1L]],
    repeats = nrow(x),
    cells_per_side = x$cells_per_side[[1L]],
    auc_median = stats::median(x$auc, na.rm = TRUE),
    auc_sd = stats::sd(x$auc, na.rm = TRUE),
    z_median = stats::median(x$z_one_sided_enriched, na.rm = TRUE),
    auc_cutoff = auc_cutoff,
    proportion_auc_ge_cutoff = mean(x$auc >= auc_cutoff, na.rm = TRUE),
    stringsAsFactors = FALSE
  ))
  do.call(rbind, rows)
}

compare_wrs_repeats <- function(first, second, top_n, auc_cutoff) {
  first <- first[order(first$pval_one_sided_enriched, -first$auc), , drop = FALSE]
  second <- second[order(second$pval_one_sided_enriched, -second$auc), , drop = FALSE]
  top_first <- utils::head(as.character(first$feature), top_n)
  top_second <- utils::head(as.character(second$feature), top_n)
  selected_first <- as.character(first$feature[first$auc >= auc_cutoff])
  selected_second <- as.character(second$feature[second$auc >= auc_cutoff])
  common <- intersect(first$feature, second$feature)
  first_common <- first[match(common, first$feature), , drop = FALSE]
  second_common <- second[match(common, second$feature), , drop = FALSE]
  top_union <- union(top_first, top_second)
  selected_union <- union(selected_first, selected_second)
  data.frame(
    repeat1 = unique(first[["repeat"]])[[1L]],
    repeat2 = unique(second[["repeat"]])[[1L]],
    n_common_features = length(common),
    top_n = top_n,
    top_overlap_n = length(intersect(top_first, top_second)),
    top_jaccard = if (length(top_union)) length(intersect(top_first, top_second)) / length(top_union) else NA_real_,
    auc_cutoff = auc_cutoff,
    selected_n_repeat1 = length(selected_first),
    selected_n_repeat2 = length(selected_second),
    selected_overlap_n = length(intersect(selected_first, selected_second)),
    selected_jaccard = if (length(selected_union)) {
      length(intersect(selected_first, selected_second)) / length(selected_union)
    } else NA_real_,
    selected_overlap_min = if (min(length(selected_first), length(selected_second)) > 0L) {
      length(intersect(selected_first, selected_second)) /
        min(length(selected_first), length(selected_second))
    } else NA_real_,
    auc_spearman = if (length(common) > 1L) suppressWarnings(stats::cor(
      first_common$auc, second_common$auc, method = "spearman", use = "complete.obs"
    )) else NA_real_,
    z_spearman = if (length(common) > 1L) suppressWarnings(stats::cor(
      first_common$z_one_sided_enriched,
      second_common$z_one_sided_enriched,
      method = "spearman", use = "complete.obs"
    )) else NA_real_,
    stringsAsFactors = FALSE
  )
}

#' Summarize repeat-to-repeat WRS feature stability
#'
#' @param result Raw repeated balanced WRS table.
#' @param top_n Number of top one-sided features compared between repeats.
#' @param auc_cutoff AUC threshold defining selected features.
#' @param stability_frequency_cutoff Minimum selected-repeat proportion for a stable feature.
#' @return List containing pairwise repeat and per-feature stability tables.
#' @export
summarize_wrs_stability <- function(result, top_n = 100L, auc_cutoff = 0.55,
                                    stability_frequency_cutoff = 0.8) {
  required <- c(
    "target_group", "repeat", "feature", "auc",
    "pval_one_sided_enriched", "z_one_sided_enriched"
  )
  missing <- setdiff(required, colnames(result))
  if (length(missing)) stop("Repeated WRS table is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  pairwise_rows <- list()
  feature_rows <- list()
  for (group in sort(unique(result$target_group))) {
    group_result <- result[result$target_group == group, , drop = FALSE]
    runs <- split(group_result, group_result[["repeat"]])
    if (length(runs) >= 2L) {
      pairs <- utils::combn(names(runs), 2L, simplify = FALSE)
      pairwise <- lapply(pairs, function(pair) {
        row <- compare_wrs_repeats(runs[[pair[[1L]]]], runs[[pair[[2L]]]], top_n, auc_cutoff)
        row$group <- group
        row
      })
      pairwise_rows[[group]] <- do.call(rbind, pairwise)
    }
    by_feature <- split(group_result, group_result$feature)
    feature_rows[[group]] <- do.call(rbind, lapply(by_feature, function(x) data.frame(
      group = group,
      feature = as.character(x$feature[[1L]]),
      repeats = nrow(x),
      selected_repeats = sum(x$auc >= auc_cutoff, na.rm = TRUE),
      selected_frequency = mean(x$auc >= auc_cutoff, na.rm = TRUE),
      stability_frequency_cutoff = stability_frequency_cutoff,
      stable = mean(x$auc >= auc_cutoff, na.rm = TRUE) >= stability_frequency_cutoff,
      auc_mean = mean(x$auc, na.rm = TRUE),
      auc_sd = stats::sd(x$auc, na.rm = TRUE),
      z_mean = mean(x$z_one_sided_enriched, na.rm = TRUE),
      z_sd = stats::sd(x$z_one_sided_enriched, na.rm = TRUE),
      stringsAsFactors = FALSE
    )))
  }
  pairwise <- if (length(pairwise_rows)) do.call(rbind, pairwise_rows) else data.frame()
  feature <- if (length(feature_rows)) do.call(rbind, feature_rows) else data.frame()
  if (nrow(feature)) feature <- feature[order(feature$group, -feature$selected_frequency, -feature$auc_mean), ]
  list(pairwise = pairwise, feature = feature)
}

wrs_to_group_list <- function(result, z_column = "z_one_sided_enriched_trunc") {
  required <- c("feature", "group", z_column)
  missing <- setdiff(required, colnames(result))
  if (length(missing)) stop("Missing WRS columns: ", paste(missing, collapse = ", "), call. = FALSE)
  result$feature <- as.character(result$feature)
  result$group <- as.character(result$group)
  split_result <- split(result, result$group)
  lapply(split_result, function(x) {
    value <- as.numeric(x[[z_column]])
    names(value) <- x$feature
    value
  })
}

#' Run section 4: SNP- and peak-level Wilcoxon analyses
#'
#' @param config Parsed EPIC2 configuration.
#' @param traits Optional trait subset.
#' @return Named list of generated files.
#' @export
run_wrs_stage <- function(config, traits = names(config$gwas$traits)) {
  output_dir <- stage_output_dir(config, "04_wrs")
  atac <- load_serialized(file.path(config$output_dir, "01_atac", "atac.qc.rds"))
  metadata <- get_atac_metadata(atac)
  outputs <- list()

  for (trait in traits) {
    matrices <- list(
      snp = load_serialized(file.path(config$output_dir, "02_accessibility", paste0(trait, ".snp_by_cell.rds"))),
      peak = load_serialized(file.path(config$output_dir, "02_accessibility", paste0(trait, ".peak_by_cell.rds")))
    )
    for (resolution in names(matrices)) {
      matrix <- matrices[[resolution]]
      cell_ids <- intersect(colnames(matrix), rownames(metadata))
      matrix <- matrix[, cell_ids, drop = FALSE]
      for (grouping in config$model$groupings) {
        if (!grouping %in% colnames(metadata)) stop("ATAC metadata column not found: ", grouping, call. = FALSE)
        groups <- setNames(metadata[cell_ids, grouping], cell_ids)
        result <- run_wrs_matrix(matrix, groups, config$wrs$truncation_quantiles)
        path <- file.path(output_dir, paste(trait, resolution, grouping, "rds", sep = "."))
        save_rds_atomic(result, path)
        write_table(result, sub("rds$", "tsv", path))
        outputs[[paste(trait, resolution, grouping, sep = ".")]] <- path

        balanced <- config$wrs$balanced
        if (isTRUE(balanced$enabled)) {
          group_counts <- table(groups)
          for (cells_per_side in balanced$cells_per_side) {
            eligible <- names(group_counts)[group_counts >= cells_per_side]
            repeated <- lapply(seq_along(eligible), function(index) {
              run_balanced_wrs_group(
                matrix, groups, eligible[[index]], cells_per_side,
                balanced$repeats, config$parameters$seed + index,
                config$wrs$truncation_quantiles
              )
            })
            repeated <- repeated[!vapply(repeated, is.null, logical(1))]
            if (length(repeated)) {
              repeated <- do.call(rbind, repeated)
              summary <- summarize_balanced_wrs(repeated, balanced$auc_cutoff)
              stability <- summarize_wrs_stability(
                repeated,
                top_n = balanced$top_n,
                auc_cutoff = balanced$auc_cutoff,
                stability_frequency_cutoff = balanced$stability_frequency_cutoff
              )
              prefix <- file.path(
                output_dir,
                paste(trait, resolution, grouping, paste0("balanced", cells_per_side), sep = ".")
              )
              save_rds_atomic(repeated, paste0(prefix, ".rds"))
              write_table(summary, paste0(prefix, ".summary.tsv"))
              write_table(stability$pairwise, paste0(prefix, ".pairwise_stability.tsv"))
              write_table(stability$feature, paste0(prefix, ".feature_stability.tsv"))
            }
          }
        }
      }
    }
    message("WRS stage complete for ", trait)
  }
  outputs
}
