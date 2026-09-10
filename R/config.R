deep_merge <- function(default, override) {
  if (!is.list(default) || !is.list(override)) return(override)
  out <- default
  for (name in names(override)) {
    out[[name]] <- if (name %in% names(default)) {
      deep_merge(default[[name]], override[[name]])
    } else {
      override[[name]]
    }
  }
  out
}

EPIC2_defaults <- function() {
  list(
    project = list(name = "EPIC2 analysis", output_dir = "results"),
    parameters = list(seed = 1000L, workers = 1L),
    genome = list(
      target = "hg19",
      chains = list(
        hg19_to_hg38 = NULL,
        hg38_to_hg19 = NULL
      ),
      unmapped_action = "drop",
      multimapping_action = "drop",
      duplicate_action = "drop"
    ),
    atac = list(
      assay = "ATAC",
      genome = "hg19",
      tissue_column = "tissue",
      group_column = "tissue_cell_type",
      attach_fragments = FALSE,
      fragment_files = NULL,
      validate_fragments = TRUE,
      compute_qc_metrics = FALSE,
      recompute_embedding = FALSE,
      qc = list(
        tss_min = 2.5,
        frip_min = 0.3,
        count_quantile_min = 0.01,
        blacklist_max = 0.02,
        nucleosome_max = 4,
        min_cells_per_group = 50L,
        peak_width_min = 100L,
        peak_width_max = 1000L,
        peak_min_total_counts = 1,
        peak_min_cells = 1L,
        peak_autosomes_only = TRUE,
        peak_remove_blacklist = TRUE
      )
    ),
    gwas = list(
      columns = list(
        rsid = "rsid", chr = "chr", pos = "pos", a1 = "A1", a2 = "A2",
        maf = "MAF", p = "P", beta = NULL, se = NULL
      ),
      maf_min = 0.05,
      genotype_format = "auto",
      genotype_orientation = "auto",
      genotype_ploidy = 2L,
      genotype_missing_values = -9,
      genotype_variant_id_column = NULL,
      genome = "hg19",
      genome_column = NULL,
      reference_missing_rate_max = 0.01,
      remove_ambiguous = TRUE,
      ld_prune = TRUE,
      ld_prune_tool = "snprelate",
      ld_prune_noninteger = "error",
      ld_method = "composite",
      ld_threshold = 0.8,
      ld_window_bp = 500000L,
      ld_slide_max_n = NULL,
      ld_prune_maf_min = 0.005,
      ld_remove_monomorphic = TRUE,
      ld_start_pos = "random.f500",
      ld_num_threads = 1L,
      traits = list()
    ),
    accessibility = list(chunks = 20L),
    ld = list(
      window_bp = 5000000L,
      ridge = 0.01,
      eigen_tolerance = 1e-8,
      correlation_block_size = 256L,
      max_snps_for_chr_ld = 15000L
    ),
    wrs = list(
      truncation_quantiles = c(0.02, 0.98),
      balanced = list(
        enabled = FALSE,
        cells_per_side = c(500L, 1000L),
        repeats = 30L,
        auc_cutoff = 0.55,
        stability_frequency_cutoff = 0.8,
        top_n = 100L
      )
    ),
    model = list(
      groupings = c("tissue_cell_type", "tissue"),
      require_all_components = TRUE
    ),
    plots = list(
      formats = c("png", "pdf"),
      dpi = 300L,
      base_size = 11,
      max_umap_cells = 100000L,
      count_bin_max = 4L,
      ld_heatmap_max_features = 250L,
      wrs_histogram_bins = 50L,
      top_features_per_group = 3L,
      pvalue_cap = 12,
      sldsc_summary = NULL,
      fail_on_error = FALSE
    )
  )
}

is_scalar_number <- function(x, lower = -Inf, upper = Inf, integer = FALSE) {
  is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= lower && x <= upper && (!integer || x == round(x))
}

is_scalar_logical <- function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x)
}

is_scalar_choice <- function(x, choices) {
  is.character(x) && length(x) == 1L && !is.na(x) && x %in% choices
}

gwas_setting_errors <- function(settings, prefix = "gwas") {
  errors <- character()
  add <- function(valid, message) {
    if (!isTRUE(valid)) errors <<- c(errors, paste0(prefix, ".", message))
  }

  genotype_formats <- c("auto", "dosage", "hard_call", "gt", "allele", "gp", "gl", "pl", "carrier")
  genotype_orientations <- c("auto", "variants_by_samples", "samples_by_variants")
  add(
    is_scalar_choice(settings$genotype_format, genotype_formats),
    paste0("genotype_format must be one of: ", paste(genotype_formats, collapse = ", "))
  )
  add(
    is_scalar_choice(settings$genotype_orientation, genotype_orientations),
    paste0("genotype_orientation must be one of: ", paste(genotype_orientations, collapse = ", "))
  )
  add(
    is_scalar_number(settings$genotype_ploidy, lower = 1, integer = TRUE),
    "genotype_ploidy must be one positive integer"
  )
  add(
    is.null(settings$genotype_variant_id_column) ||
      (is.character(settings$genotype_variant_id_column) &&
        length(settings$genotype_variant_id_column) == 1L &&
        !is.na(settings$genotype_variant_id_column) && nzchar(settings$genotype_variant_id_column)),
    "genotype_variant_id_column must be null or one non-empty column name"
  )
  add(is_genome_build(settings$genome), "genome must be hg19/GRCh37 or hg38/GRCh38")
  add(
    is.null(settings$genome_column) ||
      (is.character(settings$genome_column) && length(settings$genome_column) == 1L &&
        !is.na(settings$genome_column) && nzchar(settings$genome_column)),
    "genome_column must be null or one non-empty column name"
  )
  add(is_scalar_number(settings$maf_min, 0, 0.5), "maf_min must be one number in [0, 0.5]")
  add(
    is_scalar_number(settings$reference_missing_rate_max, 0, 1),
    "reference_missing_rate_max must be one number in [0, 1]"
  )
  add(is_scalar_logical(settings$remove_ambiguous), "remove_ambiguous must be true or false")
  add(is_scalar_logical(settings$ld_prune), "ld_prune must be true or false")
  add(
    is_scalar_choice(settings$ld_prune_tool, "snprelate"),
    "ld_prune_tool currently supports only snprelate"
  )
  add(
    is_scalar_choice(settings$ld_prune_noninteger, c("error", "round")),
    "ld_prune_noninteger must be error or round"
  )
  add(
    is_scalar_choice(settings$ld_method, c("composite", "r", "dprime", "corr")),
    "ld_method must be composite, r, dprime, or corr"
  )
  add(is_scalar_number(settings$ld_threshold, 0, 1), "ld_threshold must be one number in [0, 1]")
  add(
    is_scalar_number(settings$ld_window_bp, lower = 1, integer = TRUE),
    "ld_window_bp must be one positive integer"
  )
  add(
    is.null(settings$ld_slide_max_n) ||
      is_scalar_number(settings$ld_slide_max_n, lower = 1, integer = TRUE),
    "ld_slide_max_n must be null or one positive integer"
  )
  add(
    is_scalar_number(settings$ld_prune_maf_min, 0, 0.5),
    "ld_prune_maf_min must be one number in [0, 0.5]"
  )
  add(
    is_scalar_logical(settings$ld_remove_monomorphic),
    "ld_remove_monomorphic must be true or false"
  )
  add(
    is_scalar_choice(settings$ld_start_pos, c("random.f500", "random", "first", "last")),
    "ld_start_pos must be random.f500, random, first, or last"
  )
  add(
    is_scalar_number(settings$ld_num_threads, lower = 1, integer = TRUE),
    "ld_num_threads must be one positive integer"
  )
  errors
}

#' Read a EPIC2 YAML configuration
#'
#' @param path Path to a YAML configuration file.
#' @param check_files Whether to verify configured input files immediately.
#' @return A validated configuration list with defaults filled in.
#' @export
read_EPIC2_config <- function(path, check_files = FALSE) {
  require_namespaces("yaml", "configuration parsing")
  path <- assert_file(path, "Configuration file")
  user_config <- yaml::read_yaml(path)
  config <- deep_merge(EPIC2_defaults(), user_config)
  config$.config_file <- path
  config$.config_dir <- dirname(path)
  config$output_dir <- resolve_config_path(config, config$project$output_dir)
  validate_EPIC2_config(config, check_files = check_files)
  config
}

#' Validate a EPIC2 configuration
#'
#' @param config Configuration list returned by [read_EPIC2_config()].
#' @param check_files Whether configured input files must already exist.
#' @return The configuration, invisibly.
#' @export
validate_EPIC2_config <- function(config, check_files = FALSE) {
  errors <- character()
  if (is.null(config$atac$input)) errors <- c(errors, "atac.input is required")
  if (length(config$gwas$traits) == 0L) errors <- c(errors, "gwas.traits must contain at least one trait")
  if (is.null(names(config$gwas$traits)) || any(!nzchar(names(config$gwas$traits)))) {
    errors <- c(errors, "every gwas.traits entry must have a trait name")
  }
  if (!is_genome_build(config$atac$genome)) {
    errors <- c(errors, "atac.genome must be hg19/GRCh37 or hg38/GRCh38")
  }
  if (!is_genome_build(config$genome$target)) {
    errors <- c(errors, "genome.target must be hg19/GRCh37 or hg38/GRCh38")
  }
  for (name in c("unmapped_action", "multimapping_action", "duplicate_action")) {
    if (!is_scalar_choice(config$genome[[name]], c("drop", "error"))) {
      errors <- c(errors, paste0("genome.", name, " must be drop or error"))
    }
  }
  chains <- config$genome$chains
  if (!is.list(chains)) {
    errors <- c(errors, "genome.chains must be a mapping containing hg19_to_hg38 and hg38_to_hg19")
    chains <- list()
  }
  for (key in c("hg19_to_hg38", "hg38_to_hg19")) {
    value <- chains[[key]] %||% NULL
    if (!is.null(value) && (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value))) {
      errors <- c(errors, paste0("genome.chains.", key, " must be null or one non-empty path"))
    }
  }
  logical_atac <- c(
    "attach_fragments", "validate_fragments", "compute_qc_metrics", "recompute_embedding"
  )
  for (name in logical_atac) {
    if (!is_scalar_logical(config$atac[[name]])) {
      errors <- c(errors, paste0("atac.", name, " must be true or false"))
    }
  }
  if (isTRUE(config$atac$attach_fragments)) {
    fragment_files <- config$atac$fragment_files
    if (!is.list(fragment_files) || !length(fragment_files) || is.null(names(fragment_files)) ||
        any(!nzchar(names(fragment_files))) ||
        any(!vapply(fragment_files, function(path) {
          is.character(path) && length(path) == 1L && !is.na(path) && nzchar(path)
        }, logical(1)))) {
      errors <- c(
        errors,
        "atac.fragment_files must be a named tissue-to-path mapping when attach_fragments is true"
      )
    }
  }
  qc <- config$atac$qc
  qc_ranges <- list(
    tss_min = c(0, Inf), frip_min = c(0, 1), count_quantile_min = c(0, 1),
    blacklist_max = c(0, 1), nucleosome_max = c(0, Inf),
    peak_min_total_counts = c(0, Inf)
  )
  for (name in names(qc_ranges)) {
    range <- qc_ranges[[name]]
    if (!is_scalar_number(qc[[name]], range[[1L]], range[[2L]])) {
      errors <- c(errors, paste0("atac.qc.", name, " must be one number in [", range[[1L]], ", ", range[[2L]], "]"))
    }
  }
  qc_positive_integers <- c("min_cells_per_group", "peak_width_min", "peak_width_max", "peak_min_cells")
  for (name in qc_positive_integers) {
    if (!is_scalar_number(qc[[name]], lower = 1, integer = TRUE)) {
      errors <- c(errors, paste0("atac.qc.", name, " must be one positive integer"))
    }
  }
  if (is_scalar_number(qc$peak_width_min, lower = 1) &&
      is_scalar_number(qc$peak_width_max, lower = 1) &&
      qc$peak_width_min > qc$peak_width_max) {
    errors <- c(errors, "atac.qc.peak_width_min cannot exceed peak_width_max")
  }
  for (name in c("peak_autosomes_only", "peak_remove_blacklist")) {
    if (!is_scalar_logical(qc[[name]])) {
      errors <- c(errors, paste0("atac.qc.", name, " must be true or false"))
    }
  }

  errors <- c(errors, gwas_setting_errors(config$gwas))
  if (!is_scalar_number(config$parameters$seed, lower = 0, integer = TRUE)) {
    errors <- c(errors, "parameters.seed must be one non-negative integer")
  }
  if (!is_scalar_number(config$parameters$workers, lower = 1, integer = TRUE)) {
    errors <- c(errors, "parameters.workers must be one positive integer")
  }
  if (!is_scalar_number(config$accessibility$chunks, lower = 1, integer = TRUE)) {
    errors <- c(errors, "accessibility.chunks must be one positive integer")
  }
  if (!is_scalar_number(config$ld$window_bp, lower = 1, integer = TRUE)) {
    errors <- c(errors, "ld.window_bp must be one positive integer")
  }
  if (!is_scalar_number(config$ld$ridge, lower = 0)) {
    errors <- c(errors, "ld.ridge must be one non-negative number")
  }
  if (!is_scalar_number(config$ld$eigen_tolerance, lower = .Machine$double.eps)) {
    errors <- c(errors, "ld.eigen_tolerance must be one positive number")
  }
  if (!is_scalar_number(config$ld$max_snps_for_chr_ld, lower = 1, integer = TRUE)) {
    errors <- c(errors, "ld.max_snps_for_chr_ld must be one positive integer")
  }
  if (!is_scalar_number(config$ld$correlation_block_size, lower = 1, integer = TRUE)) {
    errors <- c(errors, "ld.correlation_block_size must be one positive integer")
  }
  truncation_quantiles <- config$wrs$truncation_quantiles
  if (length(truncation_quantiles) != 2L ||
      !is.numeric(truncation_quantiles) ||
      any(!is.finite(truncation_quantiles)) ||
      any(truncation_quantiles < 0 | truncation_quantiles > 1) ||
      truncation_quantiles[[1L]] >= truncation_quantiles[[2L]]) {
    errors <- c(errors, "wrs.truncation_quantiles must be two increasing probabilities")
  }
  balanced <- config$wrs$balanced
  if (!is_scalar_logical(balanced$enabled)) {
    errors <- c(errors, "wrs.balanced.enabled must be true or false")
  }
  if (!is.numeric(balanced$cells_per_side) || length(balanced$cells_per_side) == 0L ||
      any(!is.finite(balanced$cells_per_side)) || any(balanced$cells_per_side < 1) ||
      any(balanced$cells_per_side != round(balanced$cells_per_side))) {
    errors <- c(errors, "wrs.balanced.cells_per_side must contain positive integers")
  }
  if (!is_scalar_number(balanced$repeats, lower = 1, integer = TRUE)) {
    errors <- c(errors, "wrs.balanced.repeats must be one positive integer")
  }
  if (!is_scalar_number(balanced$auc_cutoff, 0, 1)) {
    errors <- c(errors, "wrs.balanced.auc_cutoff must be one number in [0, 1]")
  }
  if (!is_scalar_number(balanced$stability_frequency_cutoff, 0, 1)) {
    errors <- c(errors, "wrs.balanced.stability_frequency_cutoff must be one number in [0, 1]")
  }
  if (!is_scalar_number(balanced$top_n, lower = 1, integer = TRUE)) {
    errors <- c(errors, "wrs.balanced.top_n must be one positive integer")
  }
  if (!is.character(config$model$groupings) || length(config$model$groupings) == 0L ||
      any(is.na(config$model$groupings)) || any(!nzchar(config$model$groupings))) {
    errors <- c(errors, "model.groupings must contain non-empty metadata column names")
  }
  if (!is_scalar_logical(config$model$require_all_components)) {
    errors <- c(errors, "model.require_all_components must be true or false")
  }
  plots <- config$plots
  if (!is.character(plots$formats) || length(plots$formats) == 0L ||
      any(is.na(plots$formats)) || any(!plots$formats %in% c("png", "pdf"))) {
    errors <- c(errors, "plots.formats must contain png and/or pdf")
  }
  for (name in c(
    "dpi", "max_umap_cells", "count_bin_max", "ld_heatmap_max_features",
    "wrs_histogram_bins", "top_features_per_group"
  )) {
    if (!is_scalar_number(plots[[name]], lower = 1, integer = TRUE)) {
      errors <- c(errors, paste0("plots.", name, " must be one positive integer"))
    }
  }
  if (!is_scalar_number(plots$base_size, lower = 1)) {
    errors <- c(errors, "plots.base_size must be one positive number")
  }
  if (!is_scalar_number(plots$pvalue_cap, lower = 1)) {
    errors <- c(errors, "plots.pvalue_cap must be one positive number")
  }
  if (!is.null(plots$sldsc_summary) &&
      (!is.character(plots$sldsc_summary) || length(plots$sldsc_summary) != 1L ||
        is.na(plots$sldsc_summary) || !nzchar(plots$sldsc_summary))) {
    errors <- c(errors, "plots.sldsc_summary must be null or one non-empty path")
  }
  if (!is_scalar_logical(plots$fail_on_error)) {
    errors <- c(errors, "plots.fail_on_error must be true or false")
  }

  required_trait_keys <- c("summary", "genotype")
  for (trait in names(config$gwas$traits)) {
    spec <- config$gwas$traits[[trait]]
    missing <- setdiff(required_trait_keys, names(spec))
    if (length(missing)) {
      errors <- c(errors, paste0("gwas.traits.", trait, " is missing: ", paste(missing, collapse = ", ")))
    }
    effective <- deep_merge(config$gwas, spec)
    errors <- c(errors, gwas_setting_errors(effective, paste0("gwas.traits.", trait)))
  }

  if (is_genome_build(config$genome$target) && is_genome_build(config$atac$genome)) {
    target_genome <- normalize_genome_build(config$genome$target)
    atac_genome <- normalize_genome_build(config$atac$genome)
    needed_pairs <- character()
    if (target_genome != atac_genome) {
      needed_pairs <- c(needed_pairs, chain_config_key(atac_genome, target_genome))
    }
    for (trait in names(config$gwas$traits)) {
      effective <- deep_merge(config$gwas, config$gwas$traits[[trait]])
      if (is_genome_build(effective$genome)) {
        source_genome <- normalize_genome_build(effective$genome)
        if (source_genome != target_genome) {
          needed_pairs <- c(needed_pairs, chain_config_key(source_genome, target_genome))
        }
        if (target_genome != atac_genome && source_genome != atac_genome) {
          needed_pairs <- c(needed_pairs, chain_config_key(target_genome, atac_genome))
        }
      }
    }
    needed_pairs <- unique(needed_pairs)
    for (key in needed_pairs) {
      if (is.null(chains[[key]]) || !nzchar(chains[[key]])) {
        errors <- c(errors, paste0("genome.chains.", key, " is required by the configured genome conversion"))
      }
    }
  }

  if (check_files) {
    inputs <- c(config$atac$input, unlist(lapply(config$gwas$traits, `[`, required_trait_keys)))
    if (isTRUE(config$atac$attach_fragments)) {
      inputs <- c(inputs, unlist(config$atac$fragment_files, use.names = FALSE))
    }
    chain_paths <- unlist(chains, use.names = FALSE)
    chain_paths <- chain_paths[!is.na(chain_paths) & nzchar(chain_paths)]
    inputs <- c(inputs, chain_paths)
    if (!is.null(config$plots$sldsc_summary)) inputs <- c(inputs, config$plots$sldsc_summary)
    inputs <- vapply(inputs, function(x) resolve_config_path(config, x), character(1))
    missing <- inputs[!file.exists(inputs)]
    if (length(missing)) errors <- c(errors, paste0("input file not found: ", missing))
  }

  if (length(errors)) {
    stop("Invalid EPIC2 configuration:\n- ", paste(errors, collapse = "\n- "), call. = FALSE)
  }
  invisible(config)
}
