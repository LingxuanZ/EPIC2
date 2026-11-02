publication_theme <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}

safe_neglog10_plot <- function(p, cap = 12) {
  p <- suppressWarnings(as.numeric(p))
  p <- pmax(p, 1e-300)
  pmin(-log10(p), cap)
}

evenly_spaced_indices <- function(n, maximum) {
  if (n <= maximum) return(seq_len(n))
  unique(as.integer(round(seq.int(1, n, length.out = maximum))))
}

adaptive_umap_point_size <- function(cells_per_panel) {
  cells_per_panel <- suppressWarnings(as.numeric(cells_per_panel))
  if (!is.finite(cells_per_panel) || cells_per_panel < 1) return(0.15)
  max(0.15, min(2.2, 10 / sqrt(cells_per_panel)))
}

sanitize_plot_id <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  gsub("_+", "_", x)
}

publication_figure <- function(plot, width = 10, height = 7) {
  list(plot = plot, width = width, height = height)
}

save_publication_plot <- function(plot, base_path, settings, width, height) {
  require_namespaces("ggplot2", "publication figures")
  ensure_dir(dirname(base_path))
  paths <- paste0(base_path, ".", settings$formats)
  for (index in seq_along(paths)) {
    ggplot2::ggsave(
      filename = paths[[index]], plot = plot, width = width, height = height,
      units = "in", dpi = settings$dpi, bg = "white", limitsize = FALSE
    )
  }
  paths
}

new_figure_registry <- function(config, output_dir) {
  settings <- config$plots
  rows <- list()
  add <- function(id, family, status, source = NA_character_, files = character(),
                  message = NA_character_, seconds = NA_real_) {
    rows[[length(rows) + 1L]] <<- data.frame(
      figure_id = id,
      family = family,
      status = status,
      source = paste(source, collapse = ";"),
      files = paste(files, collapse = ";"),
      message = message,
      elapsed_seconds = seconds,
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }
  run <- function(id, family, subdirectory, source, builder) {
    started <- proc.time()[["elapsed"]]
    result <- tryCatch(builder(), error = identity)
    elapsed <- proc.time()[["elapsed"]] - started
    if (inherits(result, "error")) {
      add(id, family, "error", source, message = conditionMessage(result), seconds = elapsed)
      if (isTRUE(settings$fail_on_error)) stop(result)
      warning("Figure ", id, " failed: ", conditionMessage(result), call. = FALSE)
      return(invisible(NULL))
    }
    if (is.null(result) || is.null(result$plot)) {
      add(id, family, "skipped", source, message = "Required input is unavailable", seconds = elapsed)
      return(invisible(NULL))
    }
    base_path <- file.path(output_dir, subdirectory, sanitize_plot_id(id))
    paths <- tryCatch(
      save_publication_plot(result$plot, base_path, settings, result$width, result$height),
      error = identity
    )
    if (inherits(paths, "error")) {
      add(id, family, "error", source, message = conditionMessage(paths), seconds = elapsed)
      if (isTRUE(settings$fail_on_error)) stop(paths)
      warning("Figure ", id, " could not be saved: ", conditionMessage(paths), call. = FALSE)
      return(invisible(NULL))
    }
    add(id, family, "generated", source, paths, seconds = elapsed)
    invisible(paths)
  }
  skip <- function(id, family, source, message) {
    add(id, family, "skipped", source, message = message, seconds = 0)
  }
  finish <- function() {
    manifest <- if (length(rows)) do.call(rbind, rows) else data.frame(
      figure_id = character(), family = character(), status = character(),
      source = character(), files = character(), message = character(),
      elapsed_seconds = numeric(), stringsAsFactors = FALSE
    )
    write_table(manifest, file.path(output_dir, "figure_manifest.tsv"))
    manifest
  }
  list(run = run, skip = skip, finish = finish)
}

matrix_nonzero_values <- function(matrix) {
  values <- if (inherits(matrix, "sparseMatrix")) as.numeric(Matrix::summary(matrix)$x) else as.numeric(matrix)
  values[is.finite(values) & values != 0]
}

matrix_count_distribution <- function(matrix, maximum = 4L) {
  values <- matrix_nonzero_values(matrix)
  total <- as.double(nrow(matrix)) * as.double(ncol(matrix))
  bins <- c("0", as.character(seq_len(maximum)), paste0(maximum + 1L, "+"))
  counts <- setNames(numeric(length(bins)), bins)
  counts[["0"]] <- total - length(values)
  rounded <- round(values)
  integer_like <- is.finite(values) & abs(values - rounded) < 1e-8
  for (value in seq_len(maximum)) {
    counts[[as.character(value)]] <- sum(integer_like & rounded == value)
  }
  counts[[paste0(maximum + 1L, "+")]] <- sum(
    !integer_like | (integer_like & rounded > maximum)
  )
  data.frame(count = factor(names(counts), levels = names(counts)), n = unname(counts),
             proportion = unname(counts / total), stringsAsFactors = FALSE)
}

summarize_features_by_group <- function(matrix, groups, features) {
  features <- intersect(unique(as.character(features)), rownames(matrix))
  if (!length(features)) return(data.frame())
  groups <- as.character(groups[colnames(matrix)])
  indicator <- group_indicator_matrix(groups)
  sizes <- Matrix::colSums(indicator)
  subset <- matrix[features, , drop = FALSE]
  sums <- as.matrix(subset %*% indicator)
  detected <- as.matrix((subset != 0) %*% indicator)
  means <- sweep(sums, 2L, sizes, "/")
  proportions <- 100 * sweep(detected, 2L, sizes, "/")
  data.frame(
    feature = rep(rownames(subset), times = ncol(indicator)),
    cell_group = rep(colnames(indicator), each = nrow(subset)),
    pct_accessible = as.numeric(proportions),
    avg_accessibility = as.numeric(means),
    stringsAsFactors = FALSE
  )
}

standardize_within <- function(values, groups) {
  out <- numeric(length(values))
  for (group in unique(groups)) {
    index <- which(groups == group)
    x <- as.numeric(values[index])
    standard_deviation <- stats::sd(x, na.rm = TRUE)
    out[index] <- if (is.finite(standard_deviation) && standard_deviation > 0) {
      (x - mean(x, na.rm = TRUE)) / standard_deviation
    } else 0
  }
  out[!is.finite(out)] <- 0
  out
}

read_tsv_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}
