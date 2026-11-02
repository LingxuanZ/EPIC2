one_sided_p_greater <- function(statistic, degrees_freedom) {
  stats::pt(statistic, df = degrees_freedom, lower.tail = FALSE)
}

#' Fit group-wise LD-whitened linear models
#'
#' The response and complete raw design, including the intercept, are
#' transformed by the inverse square-root covariance matrix. The fit retains
#' all transformed design columns to preserve the original-scale linear model.
#'
#' @param response Named response vector.
#' @param predictors Named list of group-specific predictor vectors.
#' @param baseline Named baseline predictor vector.
#' @param inverse_sqrt Named inverse square-root covariance matrix.
#' @param trait Trait label included in output.
#' @param predictor_name Name of the focal model term.
#' @return One coefficient row per group.
#' @export
fit_whitened_groups <- function(response, predictors, baseline, inverse_sqrt,
                                trait = NA_character_, predictor_name = "predictor") {
  validate_named_matrix(inverse_sqrt, "inverse square-root covariance")
  if (nrow(inverse_sqrt) != ncol(inverse_sqrt) ||
      !identical(rownames(inverse_sqrt), colnames(inverse_sqrt))) {
    stop("inverse_sqrt must be square with identical ordered row and column names.", call. = FALSE)
  }
  if (is.null(names(response)) || is.null(names(baseline))) {
    stop("response and baseline must be named vectors.", call. = FALSE)
  }
  if (is.null(names(predictors)) || any(!nzchar(names(predictors)))) {
    stop("predictors must be a named list of named vectors.", call. = FALSE)
  }

  rows <- lapply(names(predictors), function(group) {
    predictor <- predictors[[group]]
    if (is.null(names(predictor))) stop("Predictor vector is unnamed for group: ", group, call. = FALSE)
    for (value in list(response, predictor, baseline)) {
      if (anyNA(names(value)) || any(!nzchar(names(value))) || anyDuplicated(names(value))) {
        stop("Model vectors must have unique, non-missing feature names.", call. = FALSE)
      }
    }
    ids <- rownames(inverse_sqrt)
    y <- as.numeric(response[ids])
    x <- as.numeric(predictor[ids])
    base <- as.numeric(baseline[ids])
    if (any(!is.finite(y) | !is.finite(x) | !is.finite(base))) {
      # Subsetting an inverse square root is not generally the inverse square
      # root of the corresponding marginal covariance. Do not silently do so.
      return(data.frame(trait = trait, group = group, term = predictor_name,
                        beta = NA_real_, se = NA_real_, t = NA_real_,
                        p_two_sided = NA_real_, p = NA_real_, n = length(ids),
                        df = NA_integer_, status = "missing/nonfinite inputs; rebuild covariance for retained features",
                        stringsAsFactors = FALSE))
    }
    if (length(y) <= 3L) {
      return(data.frame(
        trait = trait, group = group, term = predictor_name, beta = NA_real_,
        se = NA_real_, t = NA_real_, p_two_sided = NA_real_, p = NA_real_,
        n = length(y), df = NA_integer_, status = "insufficient observations",
        stringsAsFactors = FALSE
      ))
    }

    whitening <- inverse_sqrt[ids, ids, drop = FALSE]
    design <- cbind(intercept = 1, predictor = x, baseline = base)
    transformed_design <- as.matrix(whitening %*% design)
    transformed_response <- as.numeric(whitening %*% y)
    # Retain the transformed intercept along with both transformed predictors.
    fit_design <- transformed_design
    qr_design <- qr(fit_design)
    if (qr_design$rank < ncol(fit_design)) {
      return(data.frame(
        trait = trait, group = group, term = predictor_name, beta = NA_real_,
        se = NA_real_, t = NA_real_, p_two_sided = NA_real_, p = NA_real_,
        n = length(y), df = length(y) - qr_design$rank, status = "rank deficient",
        stringsAsFactors = FALSE
      ))
    }

    fit <- stats::lm.fit(x = fit_design, y = transformed_response)
    degrees_freedom <- length(y) - fit$rank
    residual_variance <- sum(fit$residuals^2) / degrees_freedom
    covariance <- matrix(0, ncol(fit_design), ncol(fit_design),
                          dimnames = list(colnames(fit_design), colnames(fit_design)))
    covariance[fit$qr$pivot, fit$qr$pivot] <- residual_variance * chol2inv(qr.R(fit$qr))
    standard_error <- sqrt(diag(covariance))
    statistic <- fit$coefficients / standard_error
    p_two <- 2 * stats::pt(abs(statistic), df = degrees_freedom, lower.tail = FALSE)

    data.frame(
      trait = trait,
      group = group,
      term = predictor_name,
      beta = unname(fit$coefficients[["predictor"]]),
      se = unname(standard_error[["predictor"]]),
      t = unname(statistic[["predictor"]]),
      p_two_sided = unname(p_two[["predictor"]]),
      p = one_sided_p_greater(statistic[["predictor"]], degrees_freedom),
      n = length(y),
      df = degrees_freedom,
      status = "ok",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Add within-trait multiple-testing corrections
#'
#' @param table Result data frame containing trait and P-value columns.
#' @param p_column Name of the unadjusted P-value column.
#' @return Data frame with BH and Bonferroni columns.
#' @export
adjust_pvalues_within_trait <- function(table, p_column = "p") {
  table$p_bh_within_trait <- NA_real_
  table$p_bonf_within_trait <- NA_real_
  for (trait in unique(table$trait)) {
    index <- which(table$trait == trait & is.finite(table[[p_column]]))
    table$p_bh_within_trait[index] <- stats::p.adjust(table[[p_column]][index], method = "BH")
    table$p_bonf_within_trait[index] <- stats::p.adjust(table[[p_column]][index], method = "bonferroni")
  }
  table
}

acat_pvalue <- function(pvalues) {
  pvalues <- as.numeric(pvalues)
  pvalues <- pvalues[is.finite(pvalues) & pvalues >= 0 & pvalues <= 1]
  if (length(pvalues) == 0L) return(NA_real_)
  if (any(pvalues == 0)) return(0)
  pvalues <- pmin(pvalues, 1 - 1e-15)
  terms <- ifelse(
    pvalues < 1e-16,
    1 / (pi * pvalues),
    tan((0.5 - pvalues) * pi)
  )
  statistic <- mean(terms)
  if (statistic > 1e15) return(min(1, 1 / (pi * statistic)))
  min(1, stats::pcauchy(statistic, lower.tail = FALSE))
}

#' Combine SNP/peak count and Wilcoxon component tests with ACAT
#'
#' @param snp_accessibility SNP accessibility coefficient table.
#' @param snp_wrs SNP Wilcoxon coefficient table.
#' @param peak_accessibility Peak accessibility coefficient table.
#' @param peak_wrs Peak Wilcoxon coefficient table.
#' @param require_all Require all four valid component P-values.
#' @return Trait-by-group component and combined P-value table.
#' @export
combine_component_tests <- function(snp_accessibility, snp_wrs,
                                    peak_accessibility, peak_wrs,
                                    require_all = FALSE) {
  prepare <- function(x, name) {
    required <- c("trait", "group", "p")
    missing <- setdiff(required, colnames(x))
    if (length(missing)) stop(name, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
    x <- x[, required, drop = FALSE]
    if (anyDuplicated(x[c("trait", "group")])) stop(name, " has duplicated trait/group rows.", call. = FALSE)
    colnames(x)[[3L]] <- name
    x
  }
  tables <- list(
    prepare(snp_accessibility, "p_snp_accessibility"),
    prepare(snp_wrs, "p_snp_wrs"),
    prepare(peak_accessibility, "p_peak_accessibility"),
    prepare(peak_wrs, "p_peak_wrs")
  )
  result <- Reduce(
    function(x, y) merge(x, y, by = c("trait", "group"), all = TRUE, sort = FALSE),
    tables
  )
  p_columns <- c("p_snp_accessibility", "p_snp_wrs", "p_peak_accessibility", "p_peak_wrs")
  p_matrix <- as.matrix(result[, p_columns, drop = FALSE])
  valid <- is.finite(p_matrix) & p_matrix >= 0 & p_matrix <= 1
  result$n_tests_combined <- rowSums(valid)
  p_matrix[!valid] <- NA_real_
  result$p_cauchy <- apply(p_matrix, 1L, function(row) {
    if (require_all && sum(is.finite(row)) != 4L) return(NA_real_)
    acat_pvalue(row)
  })
  result <- adjust_pvalues_within_trait(result, "p_cauchy")
  names(result)[names(result) == "p_bh_within_trait"] <- "p_cauchy_bh_within_trait"
  names(result)[names(result) == "p_bonf_within_trait"] <- "p_cauchy_bonf_within_trait"
  result[order(result$trait, result$p_cauchy), , drop = FALSE]
}

add_component_adjustment <- function(table) {
  pieces <- split(table, table$component)
  pieces <- lapply(pieces, adjust_pvalues_within_trait)
  do.call(rbind, pieces)
}

#' LD-adjusted models and Cauchy combination
#'
#' @param config Parsed EPIC2 configuration.
#' @param traits Optional trait subset.
#' @return Named list of output tables by grouping.
#' @export
run_model_stage <- function(config, traits = names(config$gwas$traits)) {
  config <- normalize_common_genome_config(config)
  lock <- acquire_stage_lock(config, "model")
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  output_dir <- stage_output_dir(config, "05_models")
  traits <- resolve_gwas_names(config, traits = traits)
  for (trait in traits) {
    validate_ld_completion(config, trait)
    validate_stage_receipt(file.path(config$output_dir, "04_wrs", paste0(trait, ".wrs_receipt.rds")),
                           wrs_stage_settings(config))
  }
  atac <- load_serialized(file.path(config$output_dir, "01_atac", "atac.qc.rds"))
  check_atac_build(atac, config$atac$assay, config$genome$build)
  metadata <- get_atac_metadata(atac)
  cell_scale <- compute_cell_scale(atac, config$atac$assay)
  all_results <- setNames(vector("list", length(config$model$groupings)), config$model$groupings)

  for (grouping in config$model$groupings) {
    if (!grouping %in% colnames(metadata)) stop("ATAC metadata column not found: ", grouping, call. = FALSE)
    component_tables <- list()
    remaining <- setdiff(names(config$gwas$traits), traits)
    old_path <- file.path(output_dir, paste0(grouping, ".component_results.rds"))
    if (file.exists(old_path)) {
      previous_traits <- intersect(unique(readRDS(old_path)$trait), remaining)
      missing_records <- previous_traits[!file.exists(file.path(output_dir,
        paste(previous_traits, grouping, "components.rds", sep = ".")))]
      if (length(missing_records)) stop("Existing results lack per-trait completion records: ",
                                        paste(missing_records, collapse = ", "),
                                        ". Rerun all traits in a new output directory.", call. = FALSE)
    }
    for (trait in remaining) {
      path <- file.path(output_dir, paste(trait, grouping, "components.rds", sep = "."))
      if (!file.exists(path)) next
      validate_ld_completion(config, trait)
      validate_stage_receipt(file.path(config$output_dir, "04_wrs", paste0(trait, ".wrs_receipt.rds")),
                             wrs_stage_settings(config))
      validate_stage_receipt(paste0(path, ".receipt.rds"), list(grouping = grouping, assay = config$atac$assay))
      component_tables[[trait]] <- readRDS(path)
    }
    for (trait in traits) {
      gwas <- load_serialized(file.path(config$output_dir, "02_gwas", paste0(trait, ".gwas.rds")))
      check_prepared_gwas_build(gwas, config$genome$build)
      snp_matrix <- load_serialized(file.path(config$output_dir, "02_accessibility", paste0(trait, ".snp_by_cell.rds")))
      peak_matrix <- load_serialized(file.path(config$output_dir, "02_accessibility", paste0(trait, ".peak_by_cell.rds")))
      check_accessibility_matrix_build(snp_matrix, config$genome$build, "SNP model input")
      check_accessibility_matrix_build(peak_matrix, config$genome$build, "peak model input")
      snp_inverse <- load_serialized(file.path(config$output_dir, "03_ld", paste0(trait, ".snp_inverse_sqrt.rds")))
      peak_model <- load_serialized(file.path(config$output_dir, "03_ld", paste0(trait, ".peak_ld_model.rds")))
      snp_wrs <- load_serialized(file.path(config$output_dir, "04_wrs", paste(trait, "snp", grouping, "rds", sep = ".")))
      peak_wrs <- load_serialized(file.path(config$output_dir, "04_wrs", paste(trait, "peak", grouping, "rds", sep = ".")))

      snp_cells <- intersect(colnames(snp_matrix), rownames(metadata))
      peak_cells <- intersect(colnames(peak_matrix), rownames(metadata))
      snp_groups <- setNames(metadata[snp_cells, grouping], snp_cells)
      peak_groups <- setNames(metadata[peak_cells, grouping], peak_cells)
      snp_accessibility <- build_group_accessibility(
        snp_matrix[, snp_cells, drop = FALSE], snp_groups, cell_scale
      )
      peak_accessibility <- build_group_accessibility(
        peak_matrix[, peak_cells, drop = FALSE], peak_groups, cell_scale
      )

      snp_response <- setNames(as.numeric(gwas$Z)^2, gwas$coord_id)
      peak_response <- peak_model$response
      snp_wrs_predictors <- wrs_to_group_list(snp_wrs)
      peak_wrs_predictors <- wrs_to_group_list(peak_wrs)

      tables <- list(
        snp_accessibility = fit_whitened_groups(
          snp_response, snp_accessibility$by_group, snp_accessibility$baseline,
          snp_inverse, trait, "A_kj"
        ),
        snp_wrs = fit_whitened_groups(
          snp_response, snp_wrs_predictors, snp_accessibility$baseline,
          snp_inverse, trait, "Z_kj"
        ),
        peak_accessibility = fit_whitened_groups(
          peak_response, peak_accessibility$by_group, peak_accessibility$baseline,
          peak_model$inverse_sqrt, trait, "A_kp"
        ),
        peak_wrs = fit_whitened_groups(
          peak_response, peak_wrs_predictors, peak_accessibility$baseline,
          peak_model$inverse_sqrt, trait, "Z_kp"
        )
      )
      for (component in names(tables)) tables[[component]]$component <- component
      component_tables[[trait]] <- do.call(rbind, tables)
      path <- file.path(output_dir, paste(trait, grouping, "components.rds", sep = "."))
      save_rds_atomic(component_tables[[trait]], path)
      dependencies <- c(file.path(config$output_dir, "03_ld", paste0(trait, ".ld_receipt.rds")),
                         file.path(config$output_dir, "04_wrs", paste0(trait, ".wrs_receipt.rds")))
      write_stage_receipt(paste0(path, ".receipt.rds"), file_fingerprints(dependencies), path,
                          list(grouping = grouping, assay = config$atac$assay))
      message("Association models complete for ", trait, " / ", grouping)
    }

    components <- do.call(rbind, component_tables)
    components <- add_component_adjustment(components)
    split_components <- split(components, components$component)
    combined <- combine_component_tests(
      split_components$snp_accessibility,
      split_components$snp_wrs,
      split_components$peak_accessibility,
      split_components$peak_wrs,
      require_all = config$model$require_all_components %||% FALSE
    )
    component_path <- file.path(output_dir, paste0(grouping, ".component_results.tsv"))
    combined_path <- file.path(output_dir, paste0(grouping, ".cauchy_results.tsv"))
    write_table(components, component_path)
    write_table(combined, combined_path)
    save_rds_atomic(components, sub("tsv$", "rds", component_path))
    save_rds_atomic(combined, sub("tsv$", "rds", combined_path))
    all_results[[grouping]] <- list(components = components, combined = combined)
  }
  all_results
}
