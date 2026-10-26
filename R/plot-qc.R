get_umap_embeddings <- function(atac) {
  if (!"umap" %in% names(atac@reductions)) return(NULL)
  require_namespaces("Seurat", "UMAP publication figures")
  embedding <- Seurat::Embeddings(atac, reduction = "umap")
  if (ncol(embedding) < 2L) return(NULL)
  embedding[, seq_len(2L), drop = FALSE]
}

plot_qc_dimensions <- function(table, base_size = 11) {
  long <- rbind(
    data.frame(stage = table$stage, metric = "Cells", value = table$cells),
    data.frame(stage = table$stage, metric = "Peaks", value = table$peaks)
  )
  long$stage <- factor(long$stage, levels = unique(table$stage))
  ggplot2::ggplot(long, ggplot2::aes(stage, value, fill = stage)) +
    ggplot2::geom_col(width = 0.68, show.legend = FALSE) +
    ggplot2::geom_text(
      ggplot2::aes(label = format(value, big.mark = ",", scientific = FALSE)),
      vjust = -0.35, size = 3.4
    ) +
    ggplot2::facet_wrap(~metric, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(before = "#5B8DB8", after = "#D1604D")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.14))) +
    publication_theme(base_size) +
    ggplot2::labs(title = "ATAC-seq quality-control retention", x = NULL, y = "Count")
}

plot_group_cell_counts <- function(table, base_size = 11) {
  table$group <- factor(table$group, levels = table$group[order(table$n_cells)])
  ggplot2::ggplot(table, ggplot2::aes(n_cells, group)) +
    ggplot2::geom_col(fill = "#3B7A78", width = 0.72) +
    publication_theme(base_size) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = max(5, base_size - 3))) +
    ggplot2::labs(title = "Retained nuclei by annotation group", x = "Nuclei", y = NULL)
}

plot_atac_qc_metrics <- function(atac, settings) {
  metadata <- get_atac_metadata(atac)
  metrics <- c("nCount_ATAC", "TSS.enrichment", "FRIP", "BlacklistRatio", "nucleosome_signal")
  metrics <- metrics[metrics %in% colnames(metadata)]
  if (!length(metrics)) return(NULL)
  index <- evenly_spaced_indices(nrow(metadata), settings$max_umap_cells)
  rows <- lapply(metrics, function(metric) {
    value <- suppressWarnings(as.numeric(metadata[index, metric]))
    if (metric == "nCount_ATAC") value <- log10(pmax(value, 0) + 1)
    data.frame(metric = metric, value = value, stringsAsFactors = FALSE)
  })
  data <- do.call(rbind, rows)
  data <- data[is.finite(data$value), , drop = FALSE]
  data$metric <- factor(
    data$metric, levels = metrics,
    labels = ifelse(metrics == "nCount_ATAC", "log10(nCount_ATAC + 1)", metrics)
  )
  ggplot2::ggplot(data, ggplot2::aes(value)) +
    ggplot2::geom_histogram(bins = 60, fill = "#5B8DB8", color = "white", linewidth = 0.15) +
    ggplot2::facet_wrap(~metric, scales = "free", ncol = 2L) +
    publication_theme(settings$base_size) +
    ggplot2::labs(title = "ATAC-seq QC metrics after filtering", x = NULL, y = "Nuclei")
}

plot_peak_detection_distribution <- function(atac, assay = "ATAC", base_size = 11) {
  counts <- get_atac_counts(atac, assay)
  data <- data.frame(
    proportion = as.numeric(Matrix::rowSums(counts != 0)) / ncol(counts)
  )
  ggplot2::ggplot(data, ggplot2::aes(proportion)) +
    ggplot2::geom_histogram(
      bins = 80, fill = "#3B7A78", color = "white", linewidth = 0.15
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 1, by = 0.25),
      labels = function(x) paste0(round(100 * x), "%")
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    publication_theme(base_size) +
    ggplot2::labs(
      title = "Peak detection across retained nuclei",
      x = "Nuclei with nonzero accessibility", y = "Peaks"
    )
}

plot_atac_density_scatter <- function(atac, settings) {
  metadata <- get_atac_metadata(atac)
  if (!all(c("nCount_ATAC", "TSS.enrichment") %in% colnames(metadata))) return(NULL)
  index <- evenly_spaced_indices(nrow(metadata), settings$max_umap_cells)
  data <- data.frame(
    count = log10(pmax(as.numeric(metadata$nCount_ATAC[index]), 0) + 1),
    tss = as.numeric(metadata$TSS.enrichment[index])
  )
  data <- data[stats::complete.cases(data), , drop = FALSE]
  ggplot2::ggplot(data, ggplot2::aes(count, tss)) +
    ggplot2::geom_bin_2d(bins = 70) +
    ggplot2::scale_fill_viridis_c(option = "C", trans = "log10") +
    publication_theme(settings$base_size) +
    ggplot2::labs(
      title = "ATAC fragment depth and TSS enrichment",
      x = "log10(nCount_ATAC + 1)", y = "TSS enrichment", fill = "Nuclei"
    )
}

plot_atac_annotation_umap <- function(atac, columns, settings) {
  embedding <- get_umap_embeddings(atac)
  if (is.null(embedding)) return(NULL)
  metadata <- get_atac_metadata(atac)
  columns <- columns[columns %in% colnames(metadata)]
  if (!length(columns)) return(NULL)
  cells <- intersect(rownames(embedding), rownames(metadata))
  index <- evenly_spaced_indices(length(cells), settings$max_umap_cells)
  cells <- cells[index]
  rows <- lapply(columns, function(column) data.frame(
    UMAP1 = embedding[cells, 1L],
    UMAP2 = embedding[cells, 2L],
    annotation = column,
    value = as.character(metadata[cells, column]),
    stringsAsFactors = FALSE
  ))
  data <- do.call(rbind, rows)
  data <- data[!is.na(data$value) & nzchar(data$value), , drop = FALSE]
  values <- sort(unique(data$value))
  palette <- setNames(grDevices::hcl.colors(length(values), "Dark 3"), values)
  point_size <- adaptive_umap_point_size(max(table(data$annotation)))
  ggplot2::ggplot(data, ggplot2::aes(UMAP1, UMAP2, color = value)) +
    ggplot2::geom_point(size = point_size, alpha = 0.8, stroke = 0) +
    ggplot2::facet_wrap(~annotation, nrow = 1L) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::coord_equal() +
    publication_theme(settings$base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(), legend.key.height = grid::unit(0.32, "cm"),
      legend.text = ggplot2::element_text(size = max(5, settings$base_size - 4))
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))) +
    ggplot2::labs(title = "Fetal chromatin accessibility atlas", x = NULL, y = NULL, color = NULL)
}

plot_gwas_filtering_waterfall <- function(table, base_size = 11) {
  stages <- c(
    n_input = "Input", n_reference_matched = "Reference matched",
    n_after_reference_missing = "Missingness QC", n_after_maf = "MAF QC",
    n_after_ambiguous = "Ambiguous removed", n_after_coordinate_validation = "Coordinates validated",
    n_after_autosomal = "Autosomes", n_after_ld_prune = "LD pruned", n_final = "Final"
  )
  available <- names(stages)[names(stages) %in% colnames(table)]
  rows <- lapply(seq_len(nrow(table)), function(index) data.frame(
    trait = table$gwas_name[[index]],
    stage = factor(unname(stages[available]), levels = unname(stages[available])),
    variants = as.numeric(table[index, available]),
    stringsAsFactors = FALSE
  ))
  data <- do.call(rbind, rows)
  ggplot2::ggplot(data, ggplot2::aes(stage, variants, color = trait, group = trait)) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    publication_theme(base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
    ggplot2::labs(
      title = "GWAS variant retention across processing steps",
      x = NULL, y = "Variants (log10 scale)", color = "GWAS"
    )
}

plot_genome_validation_qc <- function(table, base_size = 11) {
  columns <- c(
    n_duplicate_coordinates = "Duplicate coordinates"
  )
  rows <- lapply(names(columns), function(column) data.frame(
    trait = table$gwas_name,
    outcome = unname(columns[[column]]),
    variants = as.numeric(table[[column]]),
    stringsAsFactors = FALSE
  ))
  data <- do.call(rbind, rows)
  ggplot2::ggplot(data, ggplot2::aes(trait, variants, fill = outcome)) +
    ggplot2::geom_col(position = "stack") +
    publication_theme(base_size) +
    ggplot2::labs(
      title = "Duplicate variants excluded during coordinate validation",
      x = "GWAS", y = "Variants", fill = "Validation outcome"
    )
}

collect_accessibility_plot_data <- function(config, traits) {
  distributions <- list()
  coverage <- list()
  for (trait in traits) {
    for (resolution in c("snp", "peak")) {
      suffix <- if (resolution == "snp") "snp_by_cell.rds" else "peak_by_cell.rds"
      path <- file.path(config$output_dir, "02_accessibility", paste0(trait, ".", suffix))
      if (!file.exists(path)) next
      matrix <- load_serialized(path)
      check_accessibility_matrix_build(matrix, config$genome$build, "accessibility plot input")
      key <- paste(trait, resolution, sep = ".")
      distribution <- matrix_count_distribution(matrix, config$plots$count_bin_max)
      distribution$trait <- trait
      distribution$resolution <- resolution
      distributions[[key]] <- distribution
      cells_per_feature <- Matrix::rowSums(matrix != 0)
      coverage[[key]] <- data.frame(
        trait = trait,
        resolution = resolution,
        n_features = nrow(matrix),
        n_cells = ncol(matrix),
        n_features_detected_ge_1pct = sum(cells_per_feature / ncol(matrix) >= 0.01),
        proportion_features_detected_ge_1pct = mean(cells_per_feature / ncol(matrix) >= 0.01),
        median_cells_per_feature = stats::median(cells_per_feature),
        median_features_per_cell = stats::median(Matrix::colSums(matrix != 0)),
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    distributions = if (length(distributions)) do.call(rbind, distributions) else data.frame(),
    coverage = if (length(coverage)) do.call(rbind, coverage) else data.frame()
  )
}

plot_accessibility_count_distributions <- function(data, base_size = 11) {
  if (!nrow(data)) return(NULL)
  data$resolution <- factor(data$resolution, levels = c("snp", "peak"), labels = c("SNP", "Peak"))
  ggplot2::ggplot(data, ggplot2::aes(count, proportion, fill = resolution)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~trait, ncol = 4L) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
    publication_theme(base_size) +
    ggplot2::labs(
      title = "Distribution of per-feature per-cell ATAC counts",
      x = "ATAC fragment count", y = "Matrix entries", fill = "Resolution"
    )
}

plot_accessibility_coverage <- function(data, base_size = 11) {
  if (!nrow(data)) return(NULL)
  data$resolution <- factor(data$resolution, levels = c("snp", "peak"), labels = c("SNP", "Peak"))
  ggplot2::ggplot(
    data, ggplot2::aes(trait, proportion_features_detected_ge_1pct, fill = resolution)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.72) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    publication_theme(base_size) +
    ggplot2::labs(
      title = "Trait-specific features detected in at least 1% of nuclei",
      x = "GWAS", y = "Features", fill = "Resolution"
    )
}

collect_signal_umap_data <- function(config, atac, traits, resolution) {
  embedding <- get_umap_embeddings(atac)
  if (is.null(embedding)) return(data.frame())
  metadata <- get_atac_metadata(atac)
  group_column <- config$atac$group_column
  if (!group_column %in% colnames(metadata)) return(data.frame())
  cells <- intersect(rownames(embedding), colnames(atac))
  cells <- cells[evenly_spaced_indices(length(cells), config$plots$max_umap_cells)]
  rows <- list()
  for (trait in traits) {
    suffix <- if (resolution == "snp") "snp_by_cell.rds" else "peak_by_cell.rds"
    path <- file.path(config$output_dir, "02_accessibility", paste0(trait, ".", suffix))
    if (!file.exists(path)) next
    matrix <- load_serialized(path)
    use <- intersect(cells, colnames(matrix))
    if (!length(use)) next
    total <- as.numeric(Matrix::colSums(matrix[, use, drop = FALSE]))
    limits <- stats::quantile(total, c(0.1, 0.9), na.rm = TRUE, names = FALSE)
    clipped <- pmin(pmax(total, limits[[1L]]), limits[[2L]])
    span <- diff(limits)
    display <- if (is.finite(span) && span > 0) (clipped - limits[[1L]]) / span else rep(0, length(total))
    rows[[trait]] <- data.frame(
      UMAP1 = embedding[use, 1L], UMAP2 = embedding[use, 2L],
      trait = trait, group = as.character(metadata[use, group_column]),
      display = display, raw_total = total, stringsAsFactors = FALSE
    )
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

plot_signal_group_violin <- function(data, resolution, grouping, base_size = 11) {
  data <- data[!is.na(data$group) & nzchar(data$group) & is.finite(data$raw_total), , drop = FALSE]
  if (!nrow(data)) return(NULL)
  data$group <- factor(data$group, levels = sort(unique(data$group)))
  data$log_total <- log1p(data$raw_total)
  label <- if (resolution == "snp") "GWAS SNP positions" else "GWAS-overlapping peaks"
  ggplot2::ggplot(data, ggplot2::aes(group, log_total, fill = group)) +
    ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.2) +
    ggplot2::facet_wrap(~trait, ncol = 2L, scales = "free_y") +
    publication_theme(base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = max(5, base_size - 4)),
      legend.position = "none"
    ) +
    ggplot2::labs(
      title = paste("Per-nucleus accessibility at", label, "by", grouping),
      subtitle = "Violin and boxplot summaries use deterministic UMAP-matched cell sampling",
      x = NULL, y = "log1p(total ATAC fragments)"
    )
}

plot_signal_umap <- function(data, resolution, base_size = 11) {
  if (!nrow(data)) return(NULL)
  data <- data[order(data$display), , drop = FALSE]
  label <- if (resolution == "snp") "GWAS SNP positions" else "GWAS-overlapping peaks"
  point_size <- adaptive_umap_point_size(max(table(data$trait)))
  ggplot2::ggplot(data, ggplot2::aes(UMAP1, UMAP2, color = display)) +
    ggplot2::geom_point(size = point_size, alpha = 0.85, stroke = 0) +
    ggplot2::facet_wrap(~trait, ncol = 4L) +
    ggplot2::scale_color_viridis_c(option = "C", limits = c(0, 1)) +
    ggplot2::coord_equal() +
    publication_theme(base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = paste("Per-nucleus accessibility at", label),
      subtitle = "Color is scaled within trait after 10th-90th percentile clipping",
      x = NULL, y = NULL, color = "Relative\naccessibility"
    )
}
