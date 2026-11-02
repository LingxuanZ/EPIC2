model_group_order <- function(table, p_column = "p_cauchy") {
  best <- tapply(as.numeric(table[[p_column]]), as.character(table$group), function(x) {
    finite <- x[is.finite(x)]
    if (length(finite)) min(finite) else 1
  })
  names(sort(best, decreasing = TRUE))
}

model_significance_label <- function(bh, bonferroni) {
  ifelse(
    is.finite(bonferroni) & bonferroni < 0.05, "**",
    ifelse(is.finite(bh) & bh < 0.05, "*", "")
  )
}

plot_cauchy_dotplot <- function(table, grouping, cap = 12, base_size = 11) {
  required <- c(
    "trait", "group", "p_cauchy", "p_cauchy_bh_within_trait",
    "p_cauchy_bonf_within_trait", "n_tests_combined"
  )
  if (length(setdiff(required, colnames(table)))) return(NULL)
  table$group <- factor(table$group, levels = model_group_order(table))
  table$neg_log10_p <- safe_neglog10_plot(table$p_cauchy, cap)
  table$significance <- model_significance_label(
    table$p_cauchy_bh_within_trait, table$p_cauchy_bonf_within_trait
  )
  ggplot2::ggplot(
    table,
    ggplot2::aes(trait, group, color = neg_log10_p, size = n_tests_combined)
  ) +
    ggplot2::geom_point(alpha = 0.92, na.rm = TRUE) +
    ggplot2::geom_text(
      ggplot2::aes(label = significance), color = "black", size = 3, na.rm = TRUE
    ) +
    ggplot2::scale_color_gradient(low = "grey88", high = "#B2182B", limits = c(0, cap)) +
    ggplot2::scale_size_continuous(range = c(1.8, 6), breaks = 1:4) +
    publication_theme(base_size) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = max(5, base_size - 3))) +
    ggplot2::labs(
      title = paste("Cauchy-combined associations by", grouping),
      subtitle = "* within-trait BH P < 0.05; ** within-trait Bonferroni P < 0.05",
      x = "GWAS", y = NULL, color = "-log10(P)", size = "Tests\ncombined"
    )
}

cauchy_component_long <- function(table) {
  columns <- c(
    p_snp_accessibility = "SNP accessibility",
    p_snp_wrs = "SNP WRS",
    p_peak_accessibility = "Peak accessibility",
    p_peak_wrs = "Peak WRS",
    p_cauchy = "Cauchy combined"
  )
  columns <- columns[names(columns) %in% colnames(table)]
  rows <- lapply(names(columns), function(column) data.frame(
    trait = table$trait,
    group = table$group,
    test = unname(columns[[column]]),
    p_nominal = as.numeric(table[[column]]),
    stringsAsFactors = FALSE
  ))
  long <- do.call(rbind, rows)
  long$p_bh <- NA_real_
  long$p_bonferroni <- NA_real_
  split_id <- interaction(long$trait, long$test, drop = TRUE)
  for (id in unique(split_id)) {
    index <- which(split_id == id & is.finite(long$p_nominal))
    long$p_bh[index] <- stats::p.adjust(long$p_nominal[index], method = "BH")
    long$p_bonferroni[index] <- stats::p.adjust(long$p_nominal[index], method = "bonferroni")
  }
  long$test <- factor(long$test, levels = unname(columns))
  long
}

plot_component_heatmap <- function(table, grouping, adjustment = "nominal", cap = 12,
                                   base_size = 11) {
  long <- cauchy_component_long(table)
  p_column <- switch(
    adjustment, nominal = "p_nominal", bh = "p_bh", bonferroni = "p_bonferroni",
    stop("Unknown model plot adjustment: ", adjustment, call. = FALSE)
  )
  long$p_plot <- long[[p_column]]
  long$value <- safe_neglog10_plot(long$p_plot, cap)
  long$significance <- ifelse(is.finite(long$p_plot) & long$p_plot < 0.05, "*", "")
  long$group <- factor(long$group, levels = model_group_order(table))
  title_label <- switch(adjustment, nominal = "Nominal", bh = "BH-adjusted", bonferroni = "Bonferroni-adjusted")
  ggplot2::ggplot(long, ggplot2::aes(test, group, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = significance), size = 2.7) +
    ggplot2::facet_wrap(~trait, nrow = 1L) +
    ggplot2::scale_fill_gradient(low = "grey96", high = "#B2182B", limits = c(0, cap)) +
    publication_theme(base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      axis.text.y = ggplot2::element_text(size = max(5, base_size - 4))
    ) +
    ggplot2::labs(
      title = paste(title_label, "component and Cauchy P-values by", grouping),
      subtitle = "* displayed P < 0.05", x = NULL, y = NULL, fill = "-log10(P)"
    )
}

plot_component_effect_heatmap <- function(table, grouping, cap = NULL, base_size = 11) {
  required <- c("trait", "group", "component", "beta", "se", "p", "p_bh_within_trait")
  if (length(setdiff(required, colnames(table)))) return(NULL)
  table$z_effect <- as.numeric(table$beta) / as.numeric(table$se)
  if (is.null(cap)) {
    finite <- abs(table$z_effect[is.finite(table$z_effect)])
    cap <- max(1, min(8, ceiling(max(c(0, finite)))))
  }
  table$z_effect <- pmin(pmax(table$z_effect, -cap), cap)
  table$significance <- ifelse(
    is.finite(table$p_bh_within_trait) & table$p_bh_within_trait < 0.05, "*", ""
  )
  best_table <- data.frame(
    group = unique(table$group), p_cauchy = 1, stringsAsFactors = FALSE
  )
  group_min <- tapply(table$p, table$group, function(x) {
    finite <- x[is.finite(x)]
    if (length(finite)) min(finite) else 1
  })
  best_table$p_cauchy <- unname(group_min[best_table$group])
  table$group <- factor(table$group, levels = model_group_order(best_table))
  ggplot2::ggplot(table, ggplot2::aes(component, group, fill = z_effect)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = significance), size = 2.7) +
    ggplot2::facet_wrap(~trait, nrow = 1L) +
    ggplot2::scale_fill_gradient2(
      low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0,
      limits = c(-cap, cap)
    ) +
    publication_theme(base_size) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = max(8, base_size - 2)),
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
    ) +
    ggplot2::labs(
      title = paste("LD-adjusted component effect statistics by", grouping),
      subtitle = paste0("Fill is beta/SE; display range -", cap, " to ", cap,
                        "; * within-trait BH P < 0.05"),
      x = NULL, y = NULL, fill = "beta / SE"
    )
}

canonical_annotation_key <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[|/]", "_", x)
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_|_$", "", gsub("_+", "_", x))
}

find_first_column <- function(table, candidates) {
  found <- candidates[candidates %in% colnames(table)]
  if (length(found)) found[[1L]] else NULL
}

standardize_sldsc_summary <- function(table) {
  if ("Category" %in% colnames(table) && any(table$Category == "L2_1", na.rm = TRUE)) {
    table <- table[table$Category == "L2_1", , drop = FALSE]
  }
  annotation_column <- find_first_column(
    table, c("annotation_tested", "annotation", "group", "annotation_label")
  )
  if (is.null(annotation_column) || !"trait" %in% colnames(table)) {
    stop("S-LDSC summary needs trait and annotation/group columns.", call. = FALSE)
  }
  annotation <- as.character(table[[annotation_column]])
  table <- table[!is.na(annotation) & nzchar(annotation) & annotation != "All_fetal_ATAC_peaks" &
                   !is.na(table$trait) & nzchar(as.character(table$trait)), , drop = FALSE]
  annotation <- as.character(table[[annotation_column]])
  z_column <- find_first_column(table, c("coef_z", "Coefficient_z-score", "coefficient_z"))
  if (is.null(z_column)) stop("S-LDSC summary is missing coefficient z-scores.", call. = FALSE)
  p_column <- find_first_column(table, c("Coefficient_P_value", "coefficient_p", "coef_p"))
  coefficient_p <- stats::pnorm(as.numeric(table[[z_column]]), lower.tail = FALSE)
  label_column <- find_first_column(table, c("annotation_label", "label"))
  display <- if (is.null(label_column)) gsub("_", " ", annotation) else as.character(table[[label_column]])
  out <- data.frame(
    trait = as.character(table$trait),
    annotation = annotation,
    annotation_key = canonical_annotation_key(annotation),
    annotation_label = display,
    tissue = sub("_.*$", "", annotation),
    coefficient_z = as.numeric(table[[z_column]]),
    coefficient_p = coefficient_p,
    stringsAsFactors = FALSE
  )
  out$coefficient_p_supplied <- if (is.null(p_column)) NA_real_ else as.numeric(table[[p_column]])
  if (anyDuplicated(out[c("trait", "annotation_key")])) {
    stop("S-LDSC annotations are duplicated or collide after label normalization.", call. = FALSE)
  }
  enrichment_column <- find_first_column(table, c("Enrichment", "enrichment"))
  enrichment_p_column <- find_first_column(table, c("Enrichment_p", "enrichment_p"))
  out$enrichment <- if (is.null(enrichment_column)) NA_real_ else as.numeric(table[[enrichment_column]])
  out$enrichment_p <- if (is.null(enrichment_p_column)) NA_real_ else as.numeric(table[[enrichment_p_column]])
  out$coefficient_p_bh <- NA_real_
  out$enrichment_p_bh <- NA_real_
  for (trait in unique(out$trait)) {
    index <- which(out$trait == trait & is.finite(out$coefficient_p))
    out$coefficient_p_bh[index] <- stats::p.adjust(out$coefficient_p[index], method = "BH")
    index <- which(out$trait == trait & is.finite(out$enrichment_p))
    out$enrichment_p_bh[index] <- stats::p.adjust(out$enrichment_p[index], method = "BH")
  }
  out
}

plot_sldsc_coefficient_heatmap <- function(table, cap = 8, base_size = 11) {
  best <- tapply(abs(table$coefficient_z), table$annotation_label, max, na.rm = TRUE)
  table$annotation_label <- factor(table$annotation_label, levels = names(sort(best)))
  table$display_z <- pmin(pmax(table$coefficient_z, -cap), cap)
  table$significance <- ifelse(table$coefficient_p_bh < 0.05, "*", "")
  ggplot2::ggplot(table, ggplot2::aes(trait, annotation_label, fill = display_z)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = significance), size = 2.8) +
    ggplot2::scale_fill_gradient2(
      low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0,
      limits = c(-cap, cap)
    ) +
    publication_theme(base_size) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = max(5, base_size - 4))) +
    ggplot2::labs(
      title = "S-LDSC coefficient z-scores", subtitle = "* within-trait BH P < 0.05",
      x = "GWAS", y = NULL, fill = "Coefficient z"
    )
}

plot_sldsc_enrichment_dotplot <- function(table, cap = 12, base_size = 11) {
  table <- table[is.finite(table$enrichment) & is.finite(table$enrichment_p), , drop = FALSE]
  if (!nrow(table)) return(NULL)
  best <- tapply(table$enrichment_p, table$annotation_label, min, na.rm = TRUE)
  table$annotation_label <- factor(table$annotation_label, levels = names(sort(best, decreasing = TRUE)))
  table$negative_log <- safe_neglog10_plot(table$enrichment_p, cap)
  table$significance <- ifelse(table$enrichment_p_bh < 0.05, "*", "")
  ggplot2::ggplot(
    table, ggplot2::aes(trait, annotation_label, color = enrichment, size = negative_log)
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::geom_text(ggplot2::aes(label = significance), color = "black", size = 2.8) +
    ggplot2::scale_color_gradient2(low = "#3B4CC0", mid = "grey92", high = "#B40426", midpoint = 1) +
    ggplot2::scale_size_continuous(range = c(1, 6)) +
    publication_theme(base_size) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = max(5, base_size - 4))) +
    ggplot2::labs(
      title = "S-LDSC enrichment of fetal ATAC annotations",
      subtitle = "* within-trait BH P < 0.05", x = "GWAS", y = NULL,
      color = "Enrichment", size = "-log10(P)"
    )
}

prepare_cauchy_sldsc_comparison <- function(cauchy, sldsc) {
  cauchy_key <- data.frame(
    trait = as.character(cauchy$trait), group = as.character(cauchy$group),
    annotation_key = canonical_annotation_key(cauchy$group),
    cauchy_p = as.numeric(cauchy$p_cauchy),
    cauchy_bh = as.numeric(cauchy$p_cauchy_bh_within_trait),
    stringsAsFactors = FALSE
  )
  merge(cauchy_key, sldsc, by = c("trait", "annotation_key"), all = FALSE, sort = FALSE)
}

plot_cauchy_sldsc_heatmap <- function(table, cap = 12, base_size = 11) {
  if (!nrow(table)) return(NULL)
  rows <- rbind(
    data.frame(
      trait = table$trait, annotation_label = table$annotation_label,
      method = "Cauchy combined", p = table$cauchy_p, bh = table$cauchy_bh
    ),
    data.frame(
      trait = table$trait, annotation_label = table$annotation_label,
      method = "S-LDSC coefficient", p = table$coefficient_p, bh = table$coefficient_p_bh
    )
  )
  rows$value <- safe_neglog10_plot(rows$p, cap)
  rows$significance <- ifelse(is.finite(rows$bh) & rows$bh < 0.05, "*", "")
  best <- tapply(rows$p, rows$annotation_label, min, na.rm = TRUE)
  rows$annotation_label <- factor(rows$annotation_label, levels = names(sort(best, decreasing = TRUE)))
  rows$method <- factor(rows$method, levels = c("Cauchy combined", "S-LDSC coefficient"))
  ggplot2::ggplot(rows, ggplot2::aes(trait, annotation_label, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = significance), size = 2.8, na.rm = TRUE) +
    ggplot2::facet_wrap(~method, nrow = 1L) +
    ggplot2::scale_fill_gradient(low = "grey96", high = "#B2182B", limits = c(0, cap)) +
    publication_theme(base_size) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = max(5, base_size - 5))) +
    ggplot2::labs(
      title = "Cauchy combination versus S-LDSC tissue-cell-type signals",
      subtitle = "Tiles show -log10(P); * within-trait BH P < 0.05",
      x = "GWAS", y = NULL, fill = "-log10(P)"
    )
}

plot_cauchy_sldsc_scatter <- function(table, cap = 12, base_size = 11) {
  if (!nrow(table)) return(NULL)
  table$cauchy_signal <- safe_neglog10_plot(table$cauchy_p, cap)
  table$sldsc_signal <- safe_neglog10_plot(table$coefficient_p, cap)
  table$classification <- ifelse(
    table$cauchy_bh < 0.05 & table$coefficient_p_bh < 0.05, "Both",
    ifelse(table$cauchy_bh < 0.05, "Cauchy only",
           ifelse(table$coefficient_p_bh < 0.05, "S-LDSC only", "Neither"))
  )
  ggplot2::ggplot(table, ggplot2::aes(cauchy_signal, sldsc_signal, color = classification)) +
    ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(alpha = 0.85, size = 1.8, na.rm = TRUE) +
    ggplot2::facet_wrap(~trait, ncol = 4L) +
    ggplot2::scale_color_manual(values = c(
      Both = "#7F0000", `Cauchy only` = "#D95F02", `S-LDSC only` = "#1B9E77", Neither = "grey72"
    )) +
    publication_theme(base_size) +
    ggplot2::labs(
      title = "Cauchy versus S-LDSC significance",
      x = "Cauchy combined -log10(P)", y = "S-LDSC coefficient -log10(P)", color = NULL
    )
}
