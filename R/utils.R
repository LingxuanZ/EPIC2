`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

require_namespaces <- function(packages, context = NULL) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    suffix <- if (is.null(context)) "" else paste0(" for ", context)
    stop(
      "Missing required package", if (length(missing) > 1L) "s" else "",
      suffix, ": ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

assert_file <- function(path, label = "Input file") {
  if (is.null(path) || length(path) != 1L || !nzchar(path) || !file.exists(path)) {
    stop(label, " not found: ", path %||% "<missing>", call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

ensure_dir <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Unable to create directory: ", path, call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

resolve_config_path <- function(config, path) {
  if (is.null(path) || !nzchar(path)) return(path)
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) return(path)
  file.path(config$.config_dir %||% getwd(), path)
}

load_serialized <- function(path, object_name = NULL) {
  path <- assert_file(path)
  extension <- tolower(tools::file_ext(path))
  if (extension == "rds") return(readRDS(path))
  if (!extension %in% c("rda", "rdata")) {
    stop("Expected an .rds, .rda, or .RData file: ", path, call. = FALSE)
  }

  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  if (!is.null(object_name)) {
    if (!exists(object_name, envir = env, inherits = FALSE)) {
      stop("Object '", object_name, "' not found in ", path, call. = FALSE)
    }
    return(get(object_name, envir = env, inherits = FALSE))
  }
  if (length(loaded) != 1L) {
    stop(
      "File contains multiple objects; specify object_name. Objects: ",
      paste(loaded, collapse = ", "), call. = FALSE
    )
  }
  get(loaded[[1L]], envir = env, inherits = FALSE)
}

load_analysis_input <- function(path, object_name = NULL, label = "Input") {
  path <- assert_file(path, label)
  lower <- tolower(path)
  if (grepl("\\.(rds|rda|rdata)$", lower)) return(load_serialized(path, object_name))
  if (!is.null(object_name)) {
    stop(label, " object_name is only valid for .rda/.RData input.", call. = FALSE)
  }
  if (grepl("\\.csv(\\.gz)?$", lower)) {
    if (requireNamespace("data.table", quietly = TRUE)) {
      return(as.data.frame(data.table::fread(path, data.table = FALSE)))
    }
    return(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (grepl("\\.(tsv|txt)(\\.gz)?$", lower)) {
    if (requireNamespace("data.table", quietly = TRUE)) {
      return(as.data.frame(data.table::fread(path, data.table = FALSE)))
    }
    return(utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  stop(
    label, " must be .rds, .rda/.RData, .csv[.gz], .tsv[.gz], or .txt[.gz]: ", path,
    call. = FALSE
  )
}

save_rds_atomic <- function(object, path, compress = TRUE) {
  ensure_dir(dirname(path))
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, compress = compress)
  if (!file.rename(tmp, path)) stop("Unable to write output: ", path, call. = FALSE)
  invisible(path)
}

write_table <- function(x, path) {
  ensure_dir(dirname(path))
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path, sep = "\t", na = "NA")
  } else {
    utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  invisible(path)
}

normalize_chr <- function(chr, prefix = TRUE) {
  x <- sub("^chr", "", as.character(chr), ignore.case = TRUE)
  if (prefix) paste0("chr", x) else x
}

chr_to_order <- function(chr) {
  x <- toupper(sub("^chr", "", as.character(chr), ignore.case = TRUE))
  out <- suppressWarnings(as.integer(x))
  out[x == "X"] <- 23L
  out[x == "Y"] <- 24L
  out[x %in% c("M", "MT")] <- 25L
  out
}

coordinate_id <- function(chr, pos, width = 1L) {
  paste0(normalize_chr(chr), "-", as.integer(pos), "-", as.integer(pos) + width - 1L)
}

find_column <- function(x, candidates, label) {
  found <- candidates[candidates %in% colnames(x)]
  if (length(found) == 0L) {
    stop(
      "Cannot find ", label, " column. Tried: ", paste(candidates, collapse = ", "),
      ". Available columns: ", paste(colnames(x), collapse = ", "), call. = FALSE
    )
  }
  found[[1L]]
}

get_beta_se_columns <- function(gwas, columns = list()) {
  list(
    beta = columns$beta %||% find_column(
      gwas, c("Beta", "BETA", "beta", "b", "Effect", "effect"), "effect"
    ),
    se = columns$se %||% find_column(
      gwas, c("SE", "se", "StdErr", "stderr", "StdErrLogOR"), "standard error"
    )
  )
}

get_atac_counts <- function(atac, assay = "ATAC") {
  require_namespaces("Seurat", "reading the ATAC count matrix")
  tryCatch(
    Seurat::GetAssayData(atac, assay = assay, layer = "counts"),
    error = function(e) Seurat::GetAssayData(atac, assay = assay, slot = "counts")
  )
}

get_atac_metadata <- function(atac) {
  metadata <- tryCatch(atac[[]], error = function(e) NULL)
  if (is.null(metadata)) stop("Unable to read Seurat metadata.", call. = FALSE)
  metadata
}

validate_named_matrix <- function(x, label = "matrix") {
  if (!is.matrix(x) && !inherits(x, "Matrix")) stop(label, " must be a matrix.", call. = FALSE)
  if (is.null(rownames(x)) || anyDuplicated(rownames(x))) {
    stop(label, " must have unique row names.", call. = FALSE)
  }
  invisible(TRUE)
}

ordered_intersection <- function(reference, ...) {
  sets <- list(...)
  keep <- Reduce(intersect, c(list(reference), sets))
  reference[reference %in% keep]
}

stage_output_dir <- function(config, stage) {
  ensure_dir(file.path(config$output_dir, stage))
}

group_indicator_matrix <- function(groups, cell_scale = NULL) {
  require_namespaces("Matrix", "group aggregation")
  groups <- as.character(groups)
  keep <- !is.na(groups) & nzchar(groups)
  levels <- sort(unique(groups[keep]))
  column <- match(groups[keep], levels)
  if (is.null(cell_scale)) cell_scale <- rep(1, length(groups))
  Matrix::sparseMatrix(
    i = which(keep), j = column, x = as.numeric(cell_scale[keep]),
    dims = c(length(groups), length(levels)), dimnames = list(NULL, levels)
  )
}

trait_config <- function(config, trait) {
  out <- config$gwas$traits[[trait]]
  if (is.null(out)) stop("Trait not present in configuration: ", trait, call. = FALSE)
  out
}

resolve_gwas_names <- function(config, traits = NULL, gwas_name = NULL) {
  if (!is.null(traits) && !is.null(gwas_name)) {
    stop("Use either traits or gwas_name, not both.", call. = FALSE)
  }
  selected <- gwas_name %||% traits %||% names(config$gwas$traits)
  selected <- unique(as.character(selected))
  if (!length(selected) || anyNA(selected) || any(!nzchar(selected))) {
    stop("At least one non-empty GWAS name is required.", call. = FALSE)
  }
  unknown <- setdiff(selected, names(config$gwas$traits))
  if (length(unknown)) stop("Unknown GWAS name(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  selected
}
