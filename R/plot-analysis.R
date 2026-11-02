resolve_block_ids <- function(block, matrix) {
  if (is.numeric(block)) {
    index <- as.integer(block)
    index <- index[index >= 1L & index <= nrow(matrix)]
    return(rownames(matrix)[index])
  }
  intersect(as.character(block), rownames(matrix))
}

largest_ld_block <- function(matrix, blocks, maximum) {
  if (is.null(blocks) || !length(blocks) || is.null(rownames(matrix))) return(NULL)
  sizes <- vapply(blocks, length, integer(1))
  block <- blocks[[which.max(sizes)]]
  ids <- resolve_block_ids(block, matrix)
  if (!length(ids)) return(NULL)
  if (length(ids) > maximum) {
    # A contiguous central window preserves local structure; selection is independent of LD values.
    ids <- ids[floor((length(ids) - maximum) / 2) + seq_len(maximum)]
  }
  matrix[ids, ids, drop = FALSE]
}

ld_heatmap_data <- function(matrix) {
  dense <- as.matrix(matrix)
  size <- nrow(dense)
  data.frame(
    row = rep(seq_len(size), times = size),
    column = rep(seq_len(size), each = size),
    value = as.numeric(dense),
    stringsAsFactors = FALSE
  )
}

plot_ld_heatmap <- function(matrix, title, base_size = 11, fill_label = "Covariance") {
  if (is.null(matrix) || !nrow(matrix)) return(NULL)
  full <- ld_heatmap_data(matrix)
  detail <- full
  detail$value[detail$row == detail$column] <- NA_real_
  panel <- function(data, off_diagonal) {
    finite <- data$value[is.finite(data$value)]
    limits <- range(c(0, finite))
    if (diff(limits) == 0) limits <- c(0, 1)
    use_sqrt <- off_diagonal && all(finite >= 0)
    breaks <- if (use_sqrt) limits[[2L]] * c(0, 1/9, 4/9, 1) else pretty(limits, n = 4)
    ggplot2::ggplot(data, ggplot2::aes(column, row, fill = value)) +
      ggplot2::geom_raster() +
      ggplot2::scale_fill_viridis_c(
        option = "C", trans = if (use_sqrt) "sqrt" else "identity", limits = limits,
        breaks = breaks, labels = function(x) formatC(x, digits = 2, format = "g"),
        na.value = "grey85",
        guide = ggplot2::guide_colourbar(barwidth = grid::unit(4.2, "cm"),
          barheight = grid::unit(0.35, "cm"), title.position = "top", title.hjust = 0.5)
      ) +
      ggplot2::scale_y_reverse() +
      ggplot2::coord_equal(expand = FALSE) +
      publication_theme(base_size) +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(), legend.position = "bottom"
      ) +
      ggplot2::labs(
        title = paste0(title, "\n", if (off_diagonal) "Off-diagonal detail" else "Full matrix"),
        subtitle = paste(nrow(matrix), "central contiguous features;",
                         if (use_sqrt) "square-root colors" else "linear colors"),
        x = "Feature order within largest block", y = NULL, fill = fill_label
      )
  }
  combine_publication_plots(list(panel(full, FALSE), panel(detail, TRUE)), ncol = 2L)
}

collect_ld_block_summary <- function(config, traits) {
  rows <- list()
  for (trait in traits) {
    snp_blocks_path <- file.path(config$output_dir, "03_ld", paste0(trait, ".snp_blocks.rds"))
    if (file.exists(snp_blocks_path)) {
      blocks <- load_serialized(snp_blocks_path)
      rows[[paste0(trait, ".snp")]] <- data.frame(
        trait = trait, resolution = "SNP", block = names(blocks) %||% seq_along(blocks),
        n_features = vapply(blocks, length, integer(1)), stringsAsFactors = FALSE
      )
    }
    peak_path <- file.path(config$output_dir, "03_ld", paste0(trait, ".peak_ld_model.rds"))
    if (file.exists(peak_path)) {
      model <- load_serialized(peak_path)
      blocks <- model$block_indices
      rows[[paste0(trait, ".peak")]] <- data.frame(
        trait = trait, resolution = "Peak", block = names(blocks) %||% seq_along(blocks),
        n_features = vapply(blocks, length, integer(1)), stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

plot_ld_block_sizes <- function(data, base_size = 11) {
  if (!nrow(data)) return(NULL)
  ggplot2::ggplot(data, ggplot2::aes(trait, n_features, fill = resolution)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.25, position = ggplot2::position_dodge(width = 0.75)) +
    ggplot2::scale_y_log10() +
    publication_theme(base_size) +
    ggplot2::labs(
      title = "Local LD block sizes", x = "GWAS", y = "Features per block (log10)",
      fill = "Resolution"
    )
}

prepare_wrs_plot_table <- function(table) {
  required <- c("feature", "group", "auc", "pval")
  if (length(setdiff(required, colnames(table)))) return(NULL)
  if (!"pval_one_sided_enriched" %in% colnames(table)) {
    table$pval_one_sided_enriched <- ifelse(table$auc >= 0.5, table$pval / 2, 1 - table$pval / 2)
  }
  if (!"z_one_sided_enriched" %in% colnames(table)) {
    table$z_one_sided_enriched <- stats::qnorm(table$pval_one_sided_enriched, lower.tail = FALSE)
  }
  if (!"z_one_sided_enriched_trunc" %in% colnames(table)) {
    finite <- table$z_one_sided_enriched[is.finite(table$z_one_sided_enriched)]
    if (!length(finite)) stop("No finite WRS scores for plot winsorization.", call. = FALSE)
    limits <- stats::quantile(finite, c(0.02, 0.98), na.rm = TRUE, names = FALSE)
    table$z_one_sided_enriched_trunc <- pmin(
      pmax(table$z_one_sided_enriched, limits[[1L]]), limits[[2L]]
    )
  }
  table$feature <- as.character(table$feature)
  table$group <- as.character(table$group)
  table
}

summarize_wrs_groups <- function(table, trait, resolution, grouping) {
  table <- prepare_wrs_plot_table(table)
  if (is.null(table)) return(data.frame())
  pieces <- split(table, table$group)
  rows <- lapply(pieces, function(x) {
    p <- as.numeric(x$pval_one_sided_enriched)
    negative_log <- safe_neglog10_plot(p, cap = 300)
    z <- as.numeric(x$z_one_sided_enriched_trunc)
    data.frame(
      trait = trait, resolution = resolution, grouping = grouping,
      group = as.character(x$group[[1L]]), n_features = nrow(x),
      fraction_p_lt_0.05 = mean(p < 0.05, na.rm = TRUE),
      q95_neglog10_p = unname(stats::quantile(negative_log, 0.95, na.rm = TRUE)),
      q98_truncated_z = unname(stats::quantile(z, 0.98, na.rm = TRUE)),
      median_auc = stats::median(as.numeric(x$auc), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

collect_wrs_plot_data <- function(config, traits) {
  summary_rows <- list()
  tables <- list()
  for (trait in traits) {
    for (resolution in c("snp", "peak")) {
      for (grouping in config$model$groupings) {
        path <- file.path(
          config$output_dir, "04_wrs", paste(trait, resolution, grouping, "rds", sep = ".")
        )
        if (!file.exists(path)) next
        table <- prepare_wrs_plot_table(load_serialized(path))
        if (is.null(table)) next
        key <- paste(trait, resolution, grouping, sep = ".")
        tables[[key]] <- list(
          data = table, trait = trait, resolution = resolution,
          grouping = grouping, path = path
        )
        summary_rows[[key]] <- summarize_wrs_groups(table, trait, resolution, grouping)
      }
    }
  }
  list(
    tables = tables,
    summary = if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame()
  )
}

plot_wrs_histogram <- function(table, value, title, x_label, bins, base_size = 11) {
  table <- table[is.finite(as.numeric(table[[value]])) & !is.na(table$group), , drop = FALSE]
  if (!nrow(table)) return(NULL)
  table$plot_value <- as.numeric(table[[value]])
  ggplot2::ggplot(table, ggplot2::aes(x = plot_value)) +
    ggplot2::geom_histogram(bins = bins, fill = "#5B8DB8", color = "white", linewidth = 0.15) +
    ggplot2::facet_wrap(~group, scales = "free_y", ncol = 3L, labeller = ggplot2::label_wrap_gen(30)) +
    publication_theme(base_size) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = max(8, base_size - 2))) +
    ggplot2::labs(title = title, x = x_label, y = "Features")
}

plot_wrs_group_summary <- function(data, grouping, base_size = 11) {
  data <- data[data$grouping == grouping, , drop = FALSE]
  if (!nrow(data)) return(NULL)
  best <- tapply(data$q95_neglog10_p, data$group, max, na.rm = TRUE)
  data$group <- factor(data$group, levels = names(sort(best)))
  data$resolution <- factor(data$resolution, levels = c("snp", "peak"), labels = c("SNP", "Peak"))
  ggplot2::ggplot(
    data,
    ggplot2::aes(trait, group, size = fraction_p_lt_0.05, color = q95_neglog10_p)
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::facet_wrap(~resolution, nrow = 1L) +
    ggplot2::scale_size_continuous(
      range = c(0.5, 6), labels = function(x) paste0(round(100 * x), "%")
    ) +
    ggplot2::scale_color_viridis_c(option = "C") +
    publication_theme(base_size) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = max(5, base_size - 3))) +
    ggplot2::labs(
      title = paste("WRS enrichment summary:", grouping),
      x = "GWAS", y = NULL, size = "Features\nP < 0.05", color = "q95\n-log10(P)"
    )
}

top_wrs_features_per_group <- function(table, top_n) {
  table <- prepare_wrs_plot_table(table)
  table <- table[
    !is.na(table$group) & !is.na(table$feature) & is.finite(table$auc) &
      is.finite(table$pval_one_sided_enriched) & table$auc >= 0.5,
    , drop = FALSE
  ]
  if (!nrow(table)) return(data.frame())
  pieces <- split(table, table$group)
  selected <- lapply(pieces, function(x) {
    order_index <- order(x$pval_one_sided_enriched, -x$auc)
    x <- x[utils::head(order_index, top_n), , drop = FALSE]
    x$feature_rank <- seq_len(nrow(x))
    x
  })
  do.call(rbind, selected)
}

shorten_plot_label <- function(x, maximum = 24L) {
  x <- as.character(x)
  ifelse(nchar(x) > maximum, paste0(substr(x, 1L, maximum - 3L), "..."), x)
}

collect_top_wrs_dot_data <- function(config, atac, traits, resolution, grouping) {
  metadata <- get_atac_metadata(atac)
  if (!grouping %in% colnames(metadata)) return(data.frame())
  groups <- setNames(as.character(metadata[[grouping]]), rownames(metadata))
  rows <- list()
  for (trait in traits) {
    wrs_path <- file.path(
      config$output_dir, "04_wrs", paste(trait, resolution, grouping, "rds", sep = ".")
    )
    suffix <- if (resolution == "snp") "snp_by_cell.rds" else "peak_by_cell.rds"
    matrix_path <- file.path(config$output_dir, "02_accessibility", paste0(trait, ".", suffix))
    if (!file.exists(wrs_path) || !file.exists(matrix_path)) next
    top <- top_wrs_features_per_group(
      load_serialized(wrs_path), config$plots$top_features_per_group
    )
    if (!nrow(top)) next
    matrix <- load_serialized(matrix_path)
    cells <- intersect(colnames(matrix), names(groups))
    matrix <- matrix[, cells, drop = FALSE]
    group_use <- groups[cells]
    top <- top[top$feature %in% rownames(matrix), , drop = FALSE]
    if (!nrow(top)) next
    statistics <- summarize_features_by_group(matrix, group_use, top$feature)
    top$selection_id <- paste(
      trait, top$group, top$feature_rank, top$feature, sep = "__"
    )
    expanded <- merge(
      top[, c(
        "feature", "group", "feature_rank", "selection_id", "auc",
        "pval_one_sided_enriched"
      ), drop = FALSE],
      statistics, by = "feature", all.x = TRUE, sort = FALSE
    )
    expanded$trait <- trait
    expanded$resolution <- resolution
    expanded$top_group <- expanded$group
    expanded$feature_short <- shorten_plot_label(expanded$feature)
    expanded$feature_label <- paste0(
      expanded$top_group, " | top", expanded$feature_rank, ": ", expanded$feature_short
    )
    expanded$avg_accessibility_scaled <- standardize_within(
      log1p(expanded$avg_accessibility), expanded$selection_id
    )
    rows[[trait]] <- expanded
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

plot_top_wrs_dotplot <- function(data, resolution, grouping, base_size = 11) {
  if (!nrow(data)) return(NULL)
  traits <- unique(as.character(data$trait))
  group_levels <- sort(unique(c(as.character(data$top_group), as.character(data$cell_group))))
  ordered <- unique(data[order(match(data$trait, traits), match(data$top_group, group_levels),
                               data$feature_rank), c("selection_id", "feature_label")])
  labels <- setNames(as.character(ordered$feature_label), as.character(ordered$selection_id))
  data$selection_id <- factor(data$selection_id, levels = ordered$selection_id)
  data$cell_group <- factor(data$cell_group, levels = rev(group_levels))
  resolution_label <- if (resolution == "snp") "SNPs" else "GWAS-overlapping peaks"
  ggplot2::ggplot(
    data,
    ggplot2::aes(selection_id, cell_group, size = pct_accessible, color = avg_accessibility_scaled)
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::facet_wrap(~trait, scales = "free_x", ncol = 2L) +
    ggplot2::scale_x_discrete(labels = labels) +
    ggplot2::scale_size_continuous(range = c(0.25, 4.5)) +
    ggplot2::scale_color_gradient2(low = "#3B4CC0", mid = "grey92", high = "#B40426", midpoint = 0) +
    publication_theme(base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 5),
      axis.text.y = ggplot2::element_text(size = 6),
      panel.grid.major = ggplot2::element_line(color = "grey90")
    ) +
    ggplot2::labs(
      title = paste("Top positive-direction WRS-ranked", resolution_label),
      subtitle = paste0(
        "Top ", max(data$feature_rank), " per ", grouping,
        "; color is within-feature standardized log1p mean accessibility"
      ),
      x = NULL, y = NULL, size = "Percent\naccessible", color = "Scaled mean\naccessibility"
    )
}

parse_balanced_filename <- function(path, suffix) {
  pattern <- paste0("^(.+)\\.(snp|peak)\\.(.+)\\.balanced([0-9]+)\\.", suffix, "\\.tsv$")
  match <- regexec(pattern, basename(path))
  fields <- regmatches(basename(path), match)[[1L]]
  if (length(fields) != 5L) return(NULL)
  list(trait = fields[[2L]], resolution = fields[[3L]], grouping = fields[[4L]],
       cells_per_side = as.integer(fields[[5L]]))
}

collect_balanced_stability <- function(config, traits = names(config$gwas$traits)) {
  if (!isTRUE(config$wrs$balanced$enabled)) return(list(pairwise = data.frame(), feature = data.frame()))
  directory <- file.path(config$output_dir, "04_wrs")
  pair_files <- list.files(directory, pattern = "pairwise_stability\\.tsv$", full.names = TRUE)
  feature_files <- list.files(directory, pattern = "feature_stability\\.tsv$", full.names = TRUE)
  collect <- function(files, suffix) {
    rows <- lapply(files, function(path) {
      metadata <- parse_balanced_filename(path, suffix)
      if (is.null(metadata)) return(NULL)
      if (!metadata$trait %in% traits || !metadata$grouping %in% config$model$groupings ||
          !metadata$cells_per_side %in% config$wrs$balanced$cells_per_side) return(NULL)
      table <- read_tsv_if_exists(path)
      if (is.null(table)) return(NULL)
      for (name in names(metadata)) table[[name]] <- metadata[[name]]
      table
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows)) do.call(rbind, rows) else data.frame()
  }
  list(pairwise = collect(pair_files, "pairwise_stability"),
       feature = collect(feature_files, "feature_stability"))
}

plot_balanced_pairwise_stability <- function(data, base_size = 11) {
  if (!nrow(data) || !"selected_jaccard" %in% colnames(data)) return(NULL)
  data <- data[is.finite(data$selected_jaccard), , drop = FALSE]
  if (!nrow(data)) return(NULL)
  data$resolution <- factor(data$resolution, levels = c("snp", "peak"), labels = c("SNP", "Peak"))
  ggplot2::ggplot(
    data, ggplot2::aes(factor(cells_per_side), selected_jaccard, fill = resolution)
  ) +
    ggplot2::geom_boxplot(outlier.alpha = 0.18, na.rm = TRUE) +
    ggplot2::facet_wrap(~trait, ncol = 4L) +
    publication_theme(base_size) +
    ggplot2::labs(
      title = "Repeat-to-repeat stability after balanced downsampling",
      x = "Cells per side", y = "Jaccard overlap of selected features", fill = "Resolution"
    )
}

plot_balanced_feature_stability <- function(data, base_size = 11) {
  if (!nrow(data) || !"stable" %in% colnames(data)) return(NULL)
  stable <- tolower(as.character(data$stable)) %in% c("true", "t", "1")
  key <- interaction(data$trait, data$resolution, data$cells_per_side, drop = TRUE)
  count <- stats::aggregate(stable, list(key = key), sum)
  metadata <- unique(data.frame(
    key = key, trait = data$trait, resolution = data$resolution,
    cells_per_side = data$cells_per_side, stringsAsFactors = FALSE
  ))
  summary <- merge(metadata, count, by = "key", sort = FALSE)
  names(summary)[names(summary) == "x"] <- "stable_features"
  summary$resolution <- factor(summary$resolution, levels = c("snp", "peak"), labels = c("SNP", "Peak"))
  ggplot2::ggplot(
    summary, ggplot2::aes(factor(cells_per_side), stable_features, fill = resolution)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~trait, ncol = 4L, scales = "free_y") +
    publication_theme(base_size) +
    ggplot2::labs(
      title = "Stable WRS features after balanced downsampling",
      x = "Cells per side", y = "Stable features", fill = "Resolution"
    )
}
