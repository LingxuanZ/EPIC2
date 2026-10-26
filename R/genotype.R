validate_genotype_ploidy <- function(ploidy) {
  if (!is.numeric(ploidy) || length(ploidy) != 1L || !is.finite(ploidy) ||
      ploidy < 1 || abs(ploidy - round(ploidy)) > sqrt(.Machine$double.eps)) {
    stop("genotype_ploidy must be one positive integer.", call. = FALSE)
  }
  as.integer(round(ploidy))
}

genotype_missing_token <- function(x, missing_values = -9) {
  token <- toupper(trimws(as.character(x)))
  configured <- toupper(trimws(as.character(missing_values)))
  token %in% unique(c("", ".", "./.", ".|.", "NA", "NAN", "NULL", "--", configured))
}

orient_genotype_matrix <- function(x, variant_ids = NULL,
                                    orientation = c("auto", "variants_by_samples",
                                                    "samples_by_variants")) {
  orientation <- match.arg(orientation)
  if (orientation == "variants_by_samples") return(x)
  if (orientation == "samples_by_variants") return(t(x))
  if (is.null(variant_ids)) return(x)

  variant_ids <- unique(as.character(variant_ids[!is.na(variant_ids)]))
  row_hits <- if (is.null(rownames(x))) 0L else sum(rownames(x) %in% variant_ids)
  column_hits <- if (is.null(colnames(x))) 0L else sum(colnames(x) %in% variant_ids)
  if (row_hits > column_hits) return(x)
  if (column_hits > row_hits) return(t(x))
  if (row_hits > 0L) {
    stop(
      "Cannot infer genotype orientation because row and column names match the same number of GWAS variants. ",
      "Set genotype_orientation explicitly.", call. = FALSE
    )
  }
  if (is.null(rownames(x)) && nrow(x) == length(variant_ids) &&
      ncol(x) != length(variant_ids)) return(x)
  if (is.null(colnames(x)) && ncol(x) == length(variant_ids) &&
      nrow(x) != length(variant_ids)) return(t(x))
  stop(
    "Cannot infer genotype orientation: neither dimension has names matching GWAS variant IDs. ",
    "Supply variant row/column names or set genotype_orientation explicitly.", call. = FALSE
  )
}

orient_probability_array <- function(x, variant_ids = NULL,
                                      orientation = c("auto", "variants_by_samples",
                                                      "samples_by_variants")) {
  orientation <- match.arg(orientation)
  if (length(dim(x)) != 3L) stop("Genotype probability input must be a three-dimensional array.", call. = FALSE)
  if (orientation == "variants_by_samples") return(x)
  if (orientation == "samples_by_variants") return(aperm(x, c(2L, 1L, 3L)))
  probe <- matrix(0, nrow = dim(x)[1L], ncol = dim(x)[2L],
                  dimnames = dimnames(x)[1:2])
  oriented <- orient_genotype_matrix(probe, variant_ids, "auto")
  if (identical(dim(oriented), dim(probe)) &&
      identical(rownames(oriented), rownames(probe))) return(x)
  aperm(x, c(2L, 1L, 3L))
}

as_genotype_matrix <- function(x, variant_id_column = NULL) {
  if (is.data.frame(x) && !is.null(variant_id_column)) {
    if (!variant_id_column %in% colnames(x)) {
      stop("Genotype variant ID column not found: ", variant_id_column, call. = FALSE)
    }
    ids <- as.character(x[[variant_id_column]])
    x[[variant_id_column]] <- NULL
    rownames(x) <- ids
  }
  if (is.matrix(x) || inherits(x, "Matrix")) return(as.matrix(x))
  if (is.data.frame(x)) return(as.matrix(x))
  if (!is.null(dim(x)) && length(dim(x)) == 2L) {
    result <- tryCatch(as.matrix(x), error = function(e) NULL)
    if (!is.null(result)) return(result)
  }
  stop(
    "Genotype input must be matrix-like. Use a reader to extract a matrix from PLINK, VCF, or GDS files first.",
    call. = FALSE
  )
}

ensure_genotype_dimnames <- function(x, variant_ids = NULL) {
  if (is.null(rownames(x))) {
    stop(
      "Genotype variants must have row names or an explicit variant ID column. ",
      "GWAS variant_ids are orientation hints and cannot identify unnamed genotype rows. ",
      "Assign verified genotype variant IDs before input; matching dimensions are not sufficient.",
      call. = FALSE
    )
  }
  if (anyNA(rownames(x)) || any(!nzchar(rownames(x))) || anyDuplicated(rownames(x))) {
    stop("Genotype variant row names must be non-missing and unique.", call. = FALSE)
  }
  if (is.null(colnames(x))) colnames(x) <- paste0("sample_", seq_len(ncol(x)))
  if (anyNA(colnames(x)) || any(!nzchar(colnames(x))) || anyDuplicated(colnames(x))) {
    stop("Genotype sample column names must be non-missing and unique.", call. = FALSE)
  }
  x
}

detect_genotype_format <- function(x, ploidy, missing_values) {
  if (is.list(x) && !is.data.frame(x) && is.null(dim(x))) {
    if ("dosage" %in% names(x)) return("dosage")
    if (length(x) == ploidy + 1L && all(vapply(x, function(value) length(dim(value)) == 2L, logical(1)))) {
      return("gp")
    }
  }
  if (!is.null(dim(x)) && length(dim(x)) == 3L) {
    if (dim(x)[3L] != ploidy + 1L) {
      stop("A genotype probability array must have ploidy + 1 states in its third dimension.", call. = FALSE)
    }
    values <- suppressWarnings(as.numeric(x))
    values <- replace_numeric_missing(values, missing_values)
    state_matrix <- matrix(values, ncol = ploidy + 1L)
    called <- rowSums(is.na(state_matrix)) == 0L
    valid_gp <- any(called) && all(state_matrix[called, , drop = FALSE] >= 0) &&
      all(state_matrix[called, , drop = FALSE] <= 1) &&
      all(abs(rowSums(state_matrix[called, , drop = FALSE]) - 1) <= 1e-6)
    if (valid_gp) return("gp")
    stop("A probability/likelihood array is ambiguous; set genotype_format to gp, gl, or pl.", call. = FALSE)
  }
  probe <- as_genotype_matrix(x)
  if (is.logical(probe)) {
    stop(
      "Logical 0/1 genotype input is ambiguous between allele dosage and carrier status. ",
      "Set genotype_format: dosage only when it is truly allele count data.", call. = FALSE
    )
  }
  if (is.numeric(probe) || is.integer(probe)) {
    numeric_probe <- replace_numeric_missing(as.numeric(probe), missing_values)
    called <- numeric_probe[is.finite(numeric_probe)]
    if (ploidy == 2L && length(called) && all(called %in% c(0, 1))) {
      stop(
        "Numeric 0/1 genotype input is ambiguous between diploid allele dosage and carrier status. ",
        "Set genotype_format: dosage only for allele counts, or set genotype_ploidy: 1 for haploid calls.",
        call. = FALSE
      )
    }
    return("dosage")
  }

  values <- unique(as.character(probe))
  values <- values[!genotype_missing_token(values, missing_values)]
  if (!length(values)) {
    stop("Cannot infer the genotype format from an entirely missing matrix.", call. = FALSE)
  }
  numeric_values <- suppressWarnings(as.numeric(values))
  if (all(!is.na(numeric_values))) {
    if (ploidy == 2L && all(numeric_values %in% c(0, 1))) {
      stop(
        "Character 0/1 genotype input is ambiguous between diploid allele dosage and carrier status. ",
        "Set genotype_format: dosage only for allele counts, or set genotype_ploidy: 1 for haploid calls.",
        call. = FALSE
      )
    }
    return("dosage")
  }
  leading_field <- sub(":.*$", "", trimws(values))
  split_leading <- strsplit(leading_field, "[/|]")
  has_separator <- grepl("[/|]", leading_field)
  if (all(has_separator)) {
    if (all(vapply(split_leading, function(value) {
      all(value %in% c("0", "1", "."))
    }, logical(1)))) return("gt")
    return("allele")
  }
  if (any(grepl(",", values, fixed = TRUE))) {
    states <- strsplit(sub("^.*:", "", values), ",", fixed = TRUE)
    state_count <- lengths(states)
    numeric_states <- lapply(states, function(value) suppressWarnings(as.numeric(value)))
    valid_numeric <- all(vapply(numeric_states, function(value) all(!is.na(value)), logical(1)))
    valid_gp <- all(state_count == ploidy + 1L) && valid_numeric &&
      all(vapply(numeric_states, function(value) {
        all(value >= 0 & value <= 1) && abs(sum(value) - 1) <= 1e-6
      }, logical(1)))
    if (valid_gp) return("gp")
    stop("Comma-delimited genotype values are ambiguous; set genotype_format to gp, gl, or pl.", call. = FALSE)
  }

  "allele"
}

replace_numeric_missing <- function(x, missing_values) {
  if (!length(missing_values)) return(x)
  codes <- suppressWarnings(as.numeric(missing_values))
  codes <- codes[is.finite(codes)]
  for (code in unique(codes)) x[!is.na(x) & x == code] <- NA_real_
  x
}

parse_dosage_matrix <- function(x, ploidy, missing_values, hard_call = FALSE) {
  original <- x
  if (!is.numeric(x)) {
    missing <- genotype_missing_token(x, missing_values)
    values <- suppressWarnings(as.numeric(x))
    bad <- !missing & is.na(values)
    if (any(bad)) {
      example <- unique(as.character(original[bad]))[1L]
      stop("Non-numeric value in dosage genotype input: ", example, call. = FALSE)
    }
    x <- matrix(values, nrow = nrow(original), ncol = ncol(original), dimnames = dimnames(original))
    x[missing] <- NA_real_
  } else {
    storage.mode(x) <- "double"
  }
  x <- replace_numeric_missing(x, missing_values)
  if (any(!is.finite(x) & !is.na(x))) stop("Dosage genotype input contains an infinite value.", call. = FALSE)
  tolerance <- 1e-8
  if (any(x < -tolerance | x > ploidy + tolerance, na.rm = TRUE)) {
    stop("Dosages must be between 0 and genotype_ploidy after missing values are decoded.", call. = FALSE)
  }
  x[x < 0 & x >= -tolerance] <- 0
  x[x > ploidy & x <= ploidy + tolerance] <- ploidy
  if (hard_call && any(abs(x - round(x)) > tolerance, na.rm = TRUE)) {
    stop("hard_call input must contain integer allele counts.", call. = FALSE)
  }
  x
}

parse_gt_matrix <- function(x, ploidy, missing_values) {
  values <- trimws(as.character(x))
  missing <- genotype_missing_token(values, missing_values)
  gt <- sub(":.*$", "", values)
  missing <- missing | grepl("\\.", gt)
  result <- rep(NA_real_, length(values))
  called <- which(!missing)
  for (index in called) {
    alleles <- strsplit(gt[[index]], "[/|]")[[1L]]
    if (length(alleles) != ploidy || !all(alleles %in% c("0", "1"))) {
      stop(
        "GT value '", values[[index]], "' is not a biallelic ", ploidy,
        "-copy genotype. Multiallelic variants must be split before input.", call. = FALSE
      )
    }
    result[[index]] <- sum(as.integer(alleles))
  }
  matrix(result, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
}

split_allele_genotype <- function(value, ploidy) {
  value <- sub(":.*$", "", trimws(value))
  if (grepl("[/|]", value)) return(strsplit(value, "[/|]")[[1L]])
  if (nchar(value, type = "chars") == ploidy) return(strsplit(value, "", fixed = TRUE)[[1L]])
  character()
}

parse_allele_matrix <- function(x, ploidy, missing_values) {
  values <- matrix(as.character(x), nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
  result <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
  for (row in seq_len(nrow(values))) {
    parsed <- vector("list", ncol(values))
    for (column in seq_len(ncol(values))) {
      value <- values[row, column]
      if (genotype_missing_token(value, missing_values)) next
      alleles <- split_allele_genotype(value, ploidy)
      if (length(alleles) != ploidy || any(!nzchar(alleles)) || any(alleles == ".")) {
        stop(
          "Allele genotype '", value, "' at variant ", rownames(values)[[row]],
          " cannot be parsed as a ", ploidy, "-copy genotype.", call. = FALSE
        )
      }
      parsed[[column]] <- toupper(alleles)
    }
    observed <- sort(unique(unlist(parsed, use.names = FALSE)))
    if (length(observed) > 2L) {
      stop(
        "Variant ", rownames(values)[[row]], " has more than two observed alleles. ",
        "Split multiallelic variants before input.", call. = FALSE
      )
    }
    if (!length(observed)) next
    counted <- observed[[length(observed)]]
    for (column in seq_len(ncol(values))) {
      if (length(parsed[[column]])) result[row, column] <- sum(parsed[[column]] == counted)
    }
  }
  result
}

probability_list_to_array <- function(x, ploidy, variant_ids, orientation) {
  if (length(x) != ploidy + 1L ||
      !all(vapply(x, function(value) length(dim(value)) == 2L, logical(1)))) {
    stop("A genotype probability list must contain one matrix for each dosage state 0..ploidy.", call. = FALSE)
  }
  state_names <- toupper(sub("^P", "", names(x) %||% rep("", length(x))))
  if (length(state_names) && all(as.character(0:ploidy) %in% state_names)) {
    x <- x[match(as.character(0:ploidy), state_names)]
  }
  matrices <- lapply(x, function(value) {
    value <- orient_genotype_matrix(as_genotype_matrix(value), variant_ids, orientation)
    ensure_genotype_dimnames(value, variant_ids)
  })
  dimensions <- vapply(matrices, function(value) paste(dim(value), collapse = "x"), character(1))
  if (length(unique(dimensions)) != 1L ||
      !all(vapply(matrices[-1L], function(value) {
        identical(dimnames(value), dimnames(matrices[[1L]]))
      }, logical(1)))) {
    stop("All genotype probability-state matrices must have identical dimensions and names.", call. = FALSE)
  }
  array(
    unlist(lapply(matrices, as.numeric), use.names = FALSE),
    dim = c(dim(matrices[[1L]]), ploidy + 1L),
    dimnames = c(dimnames(matrices[[1L]]), list(as.character(0:ploidy)))
  )
}

probability_matrix_to_array <- function(x, ploidy, missing_values) {
  values <- as.character(x)
  result <- array(NA_real_, dim = c(nrow(x), ncol(x), ploidy + 1L),
                  dimnames = c(dimnames(x), list(as.character(0:ploidy))))
  for (index in seq_along(values)) {
    value <- trimws(values[[index]])
    if (genotype_missing_token(value, missing_values)) next
    fields <- strsplit(value, ":", fixed = TRUE)[[1L]]
    candidates <- fields[grepl(",", fields, fixed = TRUE)]
    if (length(candidates)) value <- candidates[[length(candidates)]]
    states <- strsplit(value, ",", fixed = TRUE)[[1L]]
    if (length(states) != ploidy + 1L) {
      stop("Each GP/GL/PL value must contain ploidy + 1 comma-delimited states.", call. = FALSE)
    }
    numeric_states <- suppressWarnings(as.numeric(states))
    if (anyNA(numeric_states)) stop("Non-numeric GP/GL/PL state found: ", value, call. = FALSE)
    location <- arrayInd(index, dim(x))
    result[location[[1L]], location[[2L]], ] <- numeric_states
  }
  result
}

probability_to_dosage <- function(probability, format, ploidy, missing_values) {
  if (dim(probability)[3L] != ploidy + 1L) {
    stop("Genotype probability input must contain states 0..genotype_ploidy.", call. = FALSE)
  }
  states <- lapply(seq_len(ploidy + 1L), function(index) {
    as.numeric(probability[, , index, drop = TRUE])
  })
  values <- do.call(cbind, states)
  values <- replace_numeric_missing(values, missing_values)
  called <- rowSums(is.na(values)) == 0L
  weights <- matrix(NA_real_, nrow = nrow(values), ncol = ncol(values))
  if (any(called)) {
    current <- values[called, , drop = FALSE]
    if (any(!is.finite(current))) stop(toupper(format), " input contains a non-finite state.", call. = FALSE)
    if (format == "gp") {
      if (any(current < 0)) stop("GP probabilities cannot be negative.", call. = FALSE)
      total <- rowSums(current)
      if (any(!is.finite(total) | total <= 0)) stop("Each called GP vector must have a positive sum.", call. = FALSE)
      current <- current / total
    } else if (format == "gl") {
      current <- 10^(current - apply(current, 1L, max))
      current <- current / rowSums(current)
    } else if (format == "pl") {
      current <- 10^(-(current - apply(current, 1L, min)) / 10)
      current <- current / rowSums(current)
    }
    weights[called, ] <- current
  }
  dosage <- as.numeric(weights %*% (0:ploidy))
  matrix(
    dosage, nrow = dim(probability)[1L], ncol = dim(probability)[2L],
    dimnames = dimnames(probability)[1:2]
  )
}

#' Normalize common genotype encodings to an allele-dosage matrix
#'
#' @param genotype Matrix-like genotype input, a three-dimensional probability
#'   array, or a list of probability-state matrices.
#' @param variant_ids Optional GWAS variant IDs used only to infer orientation.
#'   They never supply missing genotype identifiers; provide verified variant
#'   names on the input or an explicit `variant_id_column`.
#' @param format One of `auto`, `dosage`, `hard_call`, `gt`, `allele`, `gp`,
#'   `gl`, `pl`, or `carrier`.
#' @param ploidy Positive integer copy number.
#' @param orientation Whether variants are rows, columns, or inferred from IDs.
#' @param missing_values Numeric or character missing-value codes.
#' @param variant_id_column Optional variant ID column for data-frame input.
#' @return Numeric variant-by-sample allele-dosage matrix.
#' @export
normalize_genotype <- function(genotype, variant_ids = NULL, format = "auto", ploidy = 2,
                               orientation = "auto", missing_values = -9,
                               variant_id_column = NULL) {
  formats <- c("auto", "dosage", "hard_call", "gt", "allele", "gp", "gl", "pl", "carrier")
  format <- match.arg(tolower(format), formats)
  orientation <- match.arg(
    tolower(orientation), c("auto", "variants_by_samples", "samples_by_variants")
  )
  ploidy <- validate_genotype_ploidy(ploidy)
  if (format == "carrier") {
    stop(
      "Carrier-only 0/1 data do not contain enough information to recover allele frequency or MAF. ",
      "Provide allele counts, GT calls, or genotype probabilities.", call. = FALSE
    )
  }
  # An ID column is metadata, not an encoded sample. Remove it before format
  # detection so numeric dosage and GP tables are detected from sample values.
  if (is.data.frame(genotype) && !is.null(variant_id_column)) {
    genotype <- as_genotype_matrix(genotype, variant_id_column)
    variant_id_column <- NULL
  }
  if (format == "auto") format <- detect_genotype_format(genotype, ploidy, missing_values)

  if (format %in% c("gp", "gl", "pl")) {
    probability <- if (is.list(genotype) && !is.data.frame(genotype) && is.null(dim(genotype))) {
      probability_list_to_array(genotype, ploidy, variant_ids, orientation)
    } else if (!is.null(dim(genotype)) && length(dim(genotype)) == 3L) {
      orient_probability_array(genotype, variant_ids, orientation)
    } else {
      matrix_input <- as_genotype_matrix(genotype, variant_id_column)
      matrix_input <- orient_genotype_matrix(matrix_input, variant_ids, orientation)
      matrix_input <- ensure_genotype_dimnames(matrix_input, variant_ids)
      probability_matrix_to_array(matrix_input, ploidy, missing_values)
    }
    if (!is.numeric(probability)) {
      numeric_probability <- suppressWarnings(as.numeric(probability))
      declared_missing <- is.na(probability) |
        genotype_missing_token(probability, missing_values)
      bad <- is.na(numeric_probability) & !declared_missing
      if (any(bad)) {
        stop("Genotype probability array contains non-numeric values.", call. = FALSE)
      }
      numeric_probability[declared_missing] <- NA_real_
      probability <- array(numeric_probability, dim = dim(probability), dimnames = dimnames(probability))
    }
    dosage <- probability_to_dosage(probability, format, ploidy, missing_values)
  } else {
    if (is.list(genotype) && !is.data.frame(genotype) && "dosage" %in% names(genotype)) {
      genotype <- genotype$dosage
    }
    matrix_input <- as_genotype_matrix(genotype, variant_id_column)
    matrix_input <- orient_genotype_matrix(matrix_input, variant_ids, orientation)
    matrix_input <- ensure_genotype_dimnames(matrix_input, variant_ids)
    dosage <- switch(
      format,
      dosage = parse_dosage_matrix(matrix_input, ploidy, missing_values),
      hard_call = parse_dosage_matrix(matrix_input, ploidy, missing_values, hard_call = TRUE),
      gt = parse_gt_matrix(matrix_input, ploidy, missing_values),
      allele = parse_allele_matrix(matrix_input, ploidy, missing_values)
    )
  }
  dosage <- ensure_genotype_dimnames(dosage, variant_ids)
  attr(dosage, "genotype_format") <- format
  attr(dosage, "genotype_ploidy") <- ploidy
  dosage
}

as_snprelate_hard_calls <- function(genotype, ploidy, noninteger = c("error", "round")) {
  noninteger <- match.arg(noninteger)
  ploidy <- validate_genotype_ploidy(ploidy)
  scaled <- as.matrix(genotype) * 2 / ploidy
  tolerance <- 1e-8
  needs_rounding <- any(abs(scaled - round(scaled)) > tolerance, na.rm = TRUE)
  if (needs_rounding && noninteger == "error") {
    stop(
      "SNPRelate LD pruning requires diploid 0/1/2 hard calls. The normalized input contains ",
      "non-integer dosage. Set gwas.ld_prune_noninteger: round to opt into nearest-call ",
      "pruning, or set gwas.ld_prune: false.", call. = FALSE
    )
  }
  if (needs_rounding) {
    message("LD pruning is using nearest 0/1/2 hard calls; saved genotypes retain original dosages.")
  }
  scaled <- round(scaled)
  storage.mode(scaled) <- "integer"
  scaled
}
