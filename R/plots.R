#' Generate publication figures from completed EPIC2 stages
#'
#' Plotting is isolated from statistical computation: this stage reads frozen
#' pipeline outputs and writes PNG/PDF figures plus a machine-readable manifest.
#'
#' @param config Parsed EPIC2 configuration.
#' @param traits Optional GWAS-name subset.
#' @return The figure manifest, invisibly.
#' @export
run_plot_stage <- function(config, traits = names(config$gwas$traits)) {
  require_namespaces(c("ggplot2", "Matrix"), "publication plot stage")
  output_dir <- stage_output_dir(config, "figures")
  table_dir <- ensure_dir(file.path(output_dir, "tables"))
  registry <- new_figure_registry(config, output_dir)
  as_figure <- function(plot, width = 10, height = 7) {
    if (is.null(plot)) NULL else publication_figure(plot, width, height)
  }

  atac_path <- file.path(config$output_dir, "01_atac", "atac.qc.rds")
  atac <- if (file.exists(atac_path)) load_serialized(atac_path) else NULL
  qc_dimensions_path <- file.path(config$output_dir, "01_atac", "qc_dimensions.tsv")
  if (file.exists(qc_dimensions_path)) {
    registry$run(
      "atac_qc_dimensions", "ATAC QC", "00_atac_qc", qc_dimensions_path,
      function() as_figure(
        plot_qc_dimensions(read_tsv_if_exists(qc_dimensions_path), config$plots$base_size), 8, 5
      )
    )
  } else registry$skip("atac_qc_dimensions", "ATAC QC", qc_dimensions_path, "QC table not found")

  group_counts_path <- file.path(config$output_dir, "01_atac", "group_cell_counts.tsv")
  if (file.exists(group_counts_path)) {
    registry$run(
      "atac_group_cell_counts", "ATAC QC", "00_atac_qc", group_counts_path,
      function() {
        table <- read_tsv_if_exists(group_counts_path)
        as_figure(plot_group_cell_counts(table, config$plots$base_size), 9,
                  max(5, 0.24 * nrow(table) + 2.5))
      }
    )
  } else registry$skip("atac_group_cell_counts", "ATAC QC", group_counts_path, "Group table not found")

  if (!is.null(atac)) {
    registry$run(
      "atac_qc_metric_distributions", "ATAC QC", "00_atac_qc", atac_path,
      function() as_figure(plot_atac_qc_metrics(atac, config$plots), 10, 8)
    )
    registry$run(
      "density_scatter_nCount_TSS", "ATAC QC", "00_atac_qc", atac_path,
      function() as_figure(plot_atac_density_scatter(atac, config$plots), 7.5, 6.5)
    )
    registry$run(
      "peak_coverage_distribution", "ATAC QC", "00_atac_qc", atac_path,
      function() as_figure(
        plot_peak_detection_distribution(
          atac, config$atac$assay, config$plots$base_size
        ), 8, 6
      )
    )
    registry$run(
      "UMAP_tissue_and_tissue_cell_type.cellQC.peakQC", "Manuscript", ".", atac_path,
      function() as_figure(
        plot_atac_annotation_umap(
          atac, c(config$atac$tissue_column, config$atac$group_column), config$plots
        ),
        20, 10
      )
    )
  } else {
    registry$skip("atac_qc_metric_distributions", "ATAC QC", atac_path, "Filtered ATAC object not found")
    registry$skip("density_scatter_nCount_TSS", "ATAC QC", atac_path, "Filtered ATAC object not found")
    registry$skip(
      "UMAP_tissue_and_tissue_cell_type.cellQC.peakQC", "Manuscript", atac_path,
      "Filtered ATAC object not found"
    )
  }

  gwas_qc_path <- file.path(config$output_dir, "02_gwas", "gwas_qc_summary.tsv")
  if (file.exists(gwas_qc_path)) {
    gwas_qc <- read_tsv_if_exists(gwas_qc_path)
    gwas_qc <- gwas_qc[gwas_qc$gwas_name %in% traits, , drop = FALSE]
    registry$run(
      "gwas_filtering_waterfall", "GWAS QC", "01_gwas_qc", gwas_qc_path,
      function() as_figure(plot_gwas_filtering_waterfall(gwas_qc, config$plots$base_size), 11, 7)
    )
    registry$run(
      "gwas_liftover_exclusions", "GWAS QC", "01_gwas_qc", gwas_qc_path,
      function() as_figure(plot_liftover_qc(gwas_qc, config$plots$base_size), 9, 6)
    )
  } else {
    registry$skip("gwas_filtering_waterfall", "GWAS QC", gwas_qc_path, "GWAS QC summary not found")
    registry$skip("gwas_liftover_exclusions", "GWAS QC", gwas_qc_path, "GWAS QC summary not found")
  }

  accessibility <- collect_accessibility_plot_data(config, traits)
  if (nrow(accessibility$distributions)) {
    write_table(accessibility$distributions, file.path(table_dir, "accessibility_count_distributions.tsv"))
    write_table(accessibility$coverage, file.path(table_dir, "accessibility_coverage_summary.tsv"))
    registry$run(
      "accessibility_count_distributions", "Accessibility", "02_accessibility",
      file.path(config$output_dir, "02_accessibility"),
      function() as_figure(
        plot_accessibility_count_distributions(accessibility$distributions, config$plots$base_size),
        13, 8
      )
    )
    registry$run(
      "accessibility_coverage_summary", "Accessibility", "02_accessibility",
      file.path(table_dir, "accessibility_coverage_summary.tsv"),
      function() as_figure(plot_accessibility_coverage(accessibility$coverage, config$plots$base_size), 10, 6)
    )
  } else {
    registry$skip(
      "accessibility_count_distributions", "Accessibility",
      file.path(config$output_dir, "02_accessibility"), "No accessibility matrices found"
    )
    registry$skip(
      "accessibility_coverage_summary", "Accessibility",
      file.path(config$output_dir, "02_accessibility"), "No accessibility matrices found"
    )
  }

  for (resolution in c("snp", "peak")) {
    manuscript_id <- if (resolution == "snp") {
      "all_traits_GWAS_SNP_count_FeaturePlot"
    } else "all_traits_GWAS_peak_count_FeaturePlot"
    violin_id <- if (resolution == "snp") {
      "all_traits_GWAS_SNP_count_VlnPlot"
    } else "all_traits_GWAS_peak_count_VlnPlot"
    if (!is.null(atac)) {
      umap_data <- collect_signal_umap_data(config, atac, traits, resolution)
      registry$run(
        manuscript_id, "Manuscript", ".", file.path(config$output_dir, "02_accessibility"),
        function() as_figure(
          plot_signal_umap(umap_data, resolution, config$plots$base_size), 20, 11
        )
      )
      registry$run(
        violin_id, "Accessibility", "02_accessibility",
        c(atac_path, file.path(config$output_dir, "02_accessibility")),
        function() as_figure(
          plot_signal_group_violin(
            umap_data, resolution, config$atac$group_column, config$plots$base_size
          ), 20, max(10, 4 * ceiling(length(unique(umap_data$trait)) / 2L))
        )
      )
    } else {
      registry$skip(manuscript_id, "Manuscript", atac_path, "Filtered ATAC object not found")
      registry$skip(violin_id, "Accessibility", atac_path, "Filtered ATAC object not found")
    }
  }

  block_summary <- collect_ld_block_summary(config, traits)
  if (nrow(block_summary)) {
    write_table(block_summary, file.path(table_dir, "ld_block_summary.tsv"))
    registry$run(
      "ld_block_size_distribution", "LD diagnostics", "03_ld",
      file.path(table_dir, "ld_block_summary.tsv"),
      function() as_figure(plot_ld_block_sizes(block_summary, config$plots$base_size), 10, 6)
    )
  } else registry$skip(
    "ld_block_size_distribution", "LD diagnostics", file.path(config$output_dir, "03_ld"),
    "No LD block files found"
  )

  for (trait in traits) {
    sigma_path <- file.path(config$output_dir, "03_ld", paste0(trait, ".snp_sigma.rds"))
    blocks_path <- file.path(config$output_dir, "03_ld", paste0(trait, ".snp_blocks.rds"))
    if (file.exists(sigma_path) && file.exists(blocks_path)) {
      registry$run(
        paste0(trait, ".snp_ld_heatmap"), "LD diagnostics", "03_ld",
        c(sigma_path, blocks_path),
        function() {
          matrix <- largest_ld_block(
            load_serialized(sigma_path), load_serialized(blocks_path),
            config$plots$ld_heatmap_max_features
          )
          as_figure(plot_ld_heatmap(
            matrix, paste(trait, "SNP-level local LD covariance"), config$plots$base_size
          ), 7, 6.5)
        }
      )
    }
    peak_path <- file.path(config$output_dir, "03_ld", paste0(trait, ".peak_ld_model.rds"))
    if (file.exists(peak_path)) {
      registry$run(
        paste0(trait, ".peak_ld_heatmap"), "LD diagnostics", "03_ld", peak_path,
        function() {
          model <- load_serialized(peak_path)
          matrix <- largest_ld_block(
            model$sigma, model$block_indices, config$plots$ld_heatmap_max_features
          )
          as_figure(plot_ld_heatmap(
            matrix, paste(trait, "peak-level local covariance"), config$plots$base_size
          ), 7, 6.5)
        }
      )
    }
  }

  wrs <- collect_wrs_plot_data(config, traits)
  if (nrow(wrs$summary)) {
    write_table(wrs$summary, file.path(table_dir, "wrs_group_summary.tsv"))
    for (grouping in unique(wrs$summary$grouping)) {
      registry$run(
        paste0("wrs_group_summary.", grouping), "WRS", "04_wrs",
        file.path(table_dir, "wrs_group_summary.tsv"),
        function() as_figure(
          plot_wrs_group_summary(wrs$summary, grouping, config$plots$base_size),
          13, max(7, 0.28 * length(unique(wrs$summary$group[wrs$summary$grouping == grouping])) + 3)
        )
      )
    }
    for (key in names(wrs$tables)) {
      item <- wrs$tables[[key]]
      n_groups <- length(unique(item$data$group))
      height <- max(7, 2.2 * ceiling(n_groups / 4L))
      registry$run(
        paste0(key, ".onesided_pvalue_hist"), "WRS", "04_wrs", item$path,
        function() as_figure(plot_wrs_histogram(
          item$data, "pval_one_sided_enriched",
          paste(item$trait, item$resolution, item$grouping, "one-sided WRS P-values"),
          "One-sided enriched WRS P-value", config$plots$wrs_histogram_bins,
          config$plots$base_size
        ), 14, height)
      )
      registry$run(
        paste0(key, ".truncated_z_hist"), "WRS", "04_wrs", item$path,
        function() as_figure(plot_wrs_histogram(
          item$data, "z_one_sided_enriched_trunc",
          paste(item$trait, item$resolution, item$grouping, "truncated WRS Z-scores"),
          "One-sided enriched truncated Z", config$plots$wrs_histogram_bins,
          config$plots$base_size
        ), 14, height)
      )
    }
  } else registry$skip(
    "wrs_group_summary", "WRS", file.path(config$output_dir, "04_wrs"), "No WRS result files found"
  )

  if (!is.null(atac) && config$atac$group_column %in% config$model$groupings) {
    for (resolution in c("snp", "peak")) {
      dot_data <- collect_top_wrs_dot_data(
        config, atac, traits, resolution, config$atac$group_column
      )
      manuscript_id <- if (resolution == "snp") {
        "all_traits.top3_per_group_topSNP_wilcoxon_DotPlot"
      } else "all_traits.top3_per_group_topPeak_wilcoxon_DotPlot"
      if (nrow(dot_data)) {
        write_table(dot_data, file.path(table_dir, paste0(manuscript_id, "_data.tsv")))
      }
      trait_count <- max(1L, length(unique(dot_data$trait)))
      figure_width <- if (trait_count < 4L) 14 else 28
      figure_height <- if (trait_count < 4L) 8 else 16
      registry$run(
        manuscript_id, "Manuscript", ".", file.path(config$output_dir, "04_wrs"),
        function() as_figure(
          plot_top_wrs_dotplot(
            dot_data, resolution, config$atac$group_column, config$plots$base_size
          ),
          figure_width, figure_height
        )
      )
    }
  } else {
    registry$skip(
      "all_traits.top3_per_group_topSNP_wilcoxon_DotPlot", "Manuscript", atac_path,
      "ATAC object or configured group WRS is unavailable"
    )
    registry$skip(
      "all_traits.top3_per_group_topPeak_wilcoxon_DotPlot", "Manuscript", atac_path,
      "ATAC object or configured group WRS is unavailable"
    )
  }

  balanced <- collect_balanced_stability(config)
  if (nrow(balanced$pairwise)) {
    write_table(balanced$pairwise, file.path(table_dir, "balanced_pairwise_stability.tsv"))
    registry$run(
      "balanced_wrs_pairwise_jaccard", "WRS stability", "04_wrs",
      file.path(table_dir, "balanced_pairwise_stability.tsv"),
      function() as_figure(
        plot_balanced_pairwise_stability(balanced$pairwise, config$plots$base_size), 13, 8
      )
    )
  } else registry$skip(
    "balanced_wrs_pairwise_jaccard", "WRS stability", file.path(config$output_dir, "04_wrs"),
    "Balanced WRS pairwise files not found"
  )
  if (nrow(balanced$feature)) {
    write_table(balanced$feature, file.path(table_dir, "balanced_feature_stability.tsv"))
    registry$run(
      "balanced_wrs_stable_features", "WRS stability", "04_wrs",
      file.path(table_dir, "balanced_feature_stability.tsv"),
      function() as_figure(
        plot_balanced_feature_stability(balanced$feature, config$plots$base_size), 13, 8
      )
    )
  } else registry$skip(
    "balanced_wrs_stable_features", "WRS stability", file.path(config$output_dir, "04_wrs"),
    "Balanced WRS feature files not found"
  )

  cauchy_by_grouping <- list()
  for (grouping in config$model$groupings) {
    cauchy_path <- file.path(config$output_dir, "05_models", paste0(grouping, ".cauchy_results.tsv"))
    component_path <- file.path(
      config$output_dir, "05_models", paste0(grouping, ".component_results.tsv")
    )
    if (!file.exists(cauchy_path)) next
    cauchy <- read_tsv_if_exists(cauchy_path)
    cauchy <- cauchy[cauchy$trait %in% traits, , drop = FALSE]
    cauchy_by_grouping[[grouping]] <- cauchy
    hit <- cauchy[
      is.finite(cauchy$p_cauchy) &
        (cauchy$p_cauchy < 0.05 | cauchy$p_cauchy_bh_within_trait < 0.05),
      , drop = FALSE
    ]
    write_table(hit, file.path(table_dir, paste0(grouping, ".cauchy_nominal_or_bh_hits.tsv")))
    registry$run(
      paste0("cauchy_", grouping, "_dotplot"), "Model", "05_models", cauchy_path,
      function() as_figure(
        plot_cauchy_dotplot(
          cauchy, grouping, config$plots$pvalue_cap, config$plots$base_size
        ),
        10, max(5, 0.3 * length(unique(cauchy$group)) + 3)
      )
    )
    for (adjustment in c("nominal", "bh", "bonferroni")) {
      registry$run(
        paste0("cauchy_", grouping, "_component_heatmap_", adjustment),
        "Model", "05_models", cauchy_path,
        function() as_figure(
          plot_component_heatmap(
            cauchy, grouping, adjustment, config$plots$pvalue_cap, config$plots$base_size
          ),
          15, max(5, 0.3 * length(unique(cauchy$group)) + 3)
        )
      )
    }
    if (file.exists(component_path)) {
      components <- read_tsv_if_exists(component_path)
      components <- components[components$trait %in% traits, , drop = FALSE]
      registry$run(
        paste0("component_effects_", grouping), "Model", "05_models", component_path,
        function() as_figure(
          plot_component_effect_heatmap(components, grouping, base_size = config$plots$base_size),
          16, max(5, 0.3 * length(unique(components$group)) + 3)
        )
      )
    }
  }

  sldsc_path <- config$plots$sldsc_summary
  if (!is.null(sldsc_path)) sldsc_path <- resolve_config_path(config, sldsc_path)
  if (!is.null(sldsc_path) && file.exists(sldsc_path)) {
    sldsc <- standardize_sldsc_summary(load_analysis_input(sldsc_path, label = "S-LDSC summary"))
    write_table(sldsc, file.path(table_dir, "sldsc_standardized_summary.tsv"))
    registry$run(
      "sldsc_coefficient_zscore_heatmap", "S-LDSC", "06_sldsc", sldsc_path,
      function() as_figure(plot_sldsc_coefficient_heatmap(
        sldsc, base_size = config$plots$base_size
      ), 12, max(7, 0.28 * length(unique(sldsc$annotation_label)) + 3))
    )
    registry$run(
      "sldsc_enrichment_dotplot", "S-LDSC", "06_sldsc", sldsc_path,
      function() as_figure(plot_sldsc_enrichment_dotplot(
        sldsc, config$plots$pvalue_cap, config$plots$base_size
      ), 12, max(7, 0.28 * length(unique(sldsc$annotation_label)) + 3))
    )
    primary_grouping <- config$atac$group_column
    cauchy <- cauchy_by_grouping[[primary_grouping]]
    if (!is.null(cauchy)) {
      comparison <- prepare_cauchy_sldsc_comparison(cauchy, sldsc)
      write_table(comparison, file.path(table_dir, "cauchy_sldsc_comparison.tsv"))
      registry$run(
        "cauchy_vs_sldsc_aligned_heatmap_tissue_colored", "Manuscript", ".",
        c(file.path(config$output_dir, "05_models", paste0(primary_grouping, ".cauchy_results.tsv")),
          sldsc_path),
        function() as_figure(
          plot_cauchy_sldsc_heatmap(
            comparison, config$plots$pvalue_cap, config$plots$base_size
          ),
          17, max(8, 0.3 * length(unique(comparison$annotation_label)) + 3)
        )
      )
      registry$run(
        "cauchy_vs_sldsc_scatter", "S-LDSC", "06_sldsc",
        file.path(table_dir, "cauchy_sldsc_comparison.tsv"),
        function() as_figure(
          plot_cauchy_sldsc_scatter(
            comparison, config$plots$pvalue_cap, config$plots$base_size
          ), 13, 8
        )
      )
    } else registry$skip(
      "cauchy_vs_sldsc_aligned_heatmap_tissue_colored", "Manuscript", sldsc_path,
      "Primary grouping Cauchy table not found"
    )
  } else {
    registry$skip(
      "sldsc_coefficient_zscore_heatmap", "S-LDSC", sldsc_path %||% "<not configured>",
      "plots.sldsc_summary is not configured"
    )
    registry$skip(
      "sldsc_enrichment_dotplot", "S-LDSC", sldsc_path %||% "<not configured>",
      "plots.sldsc_summary is not configured"
    )
    registry$skip(
      "cauchy_vs_sldsc_aligned_heatmap_tissue_colored", "Manuscript",
      sldsc_path %||% "<not configured>", "plots.sldsc_summary is not configured"
    )
  }

  manifest <- registry$finish()
  generated <- sum(manifest$status == "generated")
  skipped <- sum(manifest$status == "skipped")
  failed <- sum(manifest$status == "error")
  message(
    "Plot stage complete: ", generated, " generated, ", skipped, " skipped, ", failed,
    " failed. Manifest: ", file.path(output_dir, "figure_manifest.tsv")
  )
  invisible(manifest)
}
