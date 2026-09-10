genome_build_aliases <- function() {
  c(
    hg19 = "hg19", grch37 = "hg19", b37 = "hg19", build37 = "hg19",
    hg38 = "hg38", grch38 = "hg38", b38 = "hg38", build38 = "hg38"
  )
}

#' Normalize a human genome-build name
#'
#' @param genome Genome build name such as `hg19`, `GRCh37`, `hg38`, or `GRCh38`.
#' @return Canonical UCSC build name, either `hg19` or `hg38`.
#' @export
normalize_genome_build <- function(genome) {
  if (!is.character(genome) || length(genome) != 1L || is.na(genome) || !nzchar(genome)) {
    stop("Genome build must be one non-empty value: hg19/GRCh37 or hg38/GRCh38.", call. = FALSE)
  }
  key <- tolower(gsub("[^A-Za-z0-9]", "", genome))
  value <- unname(genome_build_aliases()[key])
  if (length(value) != 1L || is.na(value)) {
    stop(
      "Unsupported genome build '", genome,
      "'. Use hg19/GRCh37 or hg38/GRCh38; the software does not guess a build from coordinates.",
      call. = FALSE
    )
  }
  value
}

is_genome_build <- function(genome) {
  !inherits(try(normalize_genome_build(genome), silent = TRUE), "try-error")
}

genome_grch_name <- function(genome) {
  if (normalize_genome_build(genome) == "hg19") "GRCh37" else "GRCh38"
}

human_genome_lengths <- function(genome) {
  genome <- normalize_genome_build(genome)
  lengths <- if (genome == "hg19") {
    c(
      249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
      159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
      115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
      59128983, 63025520, 48129895, 51304566, 155270560, 59373566, 16571
    )
  } else {
    c(
      248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
      159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
      114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
      58617616, 64444167, 46709983, 50818468, 156040895, 57227415, 16569
    )
  }
  names(lengths) <- paste0("chr", c(1:22, "X", "Y", "M"))
  lengths
}

canonical_chr <- function(chr) {
  value <- normalize_chr(chr)
  value[toupper(value) == "CHRMT"] <- "chrM"
  value
}

#' Validate human genomic coordinates against a declared build
#'
#' Coordinate bounds can detect impossible declarations, but cannot distinguish
#' hg19 from hg38 for positions valid in both builds. Input builds must therefore
#' be declared explicitly.
#'
#' @param chr Chromosome values.
#' @param start One-based inclusive starts.
#' @param end One-based inclusive ends. Defaults to `start`.
#' @param genome Declared input build.
#' @param label Input label used in errors.
#' @return A one-row coordinate QC data frame.
#' @export
validate_genome_coordinates <- function(chr, start, end = start, genome, label = "coordinates") {
  genome <- normalize_genome_build(genome)
  if (length(chr) != length(start) || length(start) != length(end)) {
    stop(label, " chromosome, start, and end vectors must have equal lengths.", call. = FALSE)
  }
  chromosome <- canonical_chr(chr)
  start <- suppressWarnings(as.integer(start))
  end <- suppressWarnings(as.integer(end))
  basic_invalid <- is.na(chromosome) | !nzchar(chromosome) |
    !is.finite(start) | !is.finite(end) | start < 1L | end < start
  lengths <- human_genome_lengths(genome)
  known <- chromosome %in% names(lengths)
  out_of_bounds <- !basic_invalid & known & end > unname(lengths[chromosome])
  if (any(basic_invalid | out_of_bounds)) {
    index <- which(basic_invalid | out_of_bounds)
    examples <- utils::head(
      paste0(chromosome[index], ":", start[index], "-", end[index]), 5L
    )
    stop(
      label, " contains ", length(index), " invalid or out-of-bounds range(s) for ",
      genome, "/", genome_grch_name(genome), ". Examples: ",
      paste(examples, collapse = ", "),
      ". This may indicate an incorrectly declared genome build.",
      call. = FALSE
    )
  }
  data.frame(
    genome = genome,
    n_ranges = length(start),
    n_standard_chromosome = sum(known),
    n_nonstandard_chromosome = sum(!known),
    n_out_of_bounds = 0L,
    stringsAsFactors = FALSE
  )
}

normalize_granges_seqnames <- function(ranges) {
  require_namespaces(c("GenomicRanges", "IRanges", "S4Vectors"), "genomic range normalization")
  result <- GenomicRanges::GRanges(
    seqnames = canonical_chr(as.character(GenomicRanges::seqnames(ranges))),
    ranges = IRanges::IRanges(
      start = GenomicRanges::start(ranges),
      end = GenomicRanges::end(ranges)
    ),
    strand = as.character(GenomicRanges::strand(ranges))
  )
  names(result) <- names(ranges)
  if (ncol(S4Vectors::mcols(ranges))) S4Vectors::mcols(result) <- S4Vectors::mcols(ranges)
  result
}

declared_range_genome <- function(ranges) {
  require_namespaces("GenomeInfoDb", "stored genome-build inspection")
  values <- unique(as.character(GenomeInfoDb::genome(ranges)))
  values <- values[!is.na(values) & nzchar(values)]
  if (!length(values)) return(NA_character_)
  recognized <- vapply(values, function(value) {
    tryCatch(normalize_genome_build(value), error = function(e) NA_character_)
  }, character(1))
  recognized <- unique(recognized[!is.na(recognized)])
  if (length(recognized) > 1L) {
    stop("Stored genomic ranges contain conflicting hg19 and hg38 genome metadata.", call. = FALSE)
  }
  recognized %||% NA_character_
}

check_range_genome_declaration <- function(ranges, declared_genome, label) {
  declared_genome <- normalize_genome_build(declared_genome)
  stored <- declared_range_genome(ranges)
  if (!is.na(stored) && stored != declared_genome) {
    stop(
      label, " is declared as ", declared_genome, " but its stored Seqinfo genome is ", stored, ".",
      call. = FALSE
    )
  }
  data.frame(
    input = label,
    declared_genome = declared_genome,
    stored_genome = stored,
    metadata_status = if (is.na(stored)) "not_available" else "matched",
    stringsAsFactors = FALSE
  )
}

check_table_genome_declaration <- function(table, declared_genome, genome_column = NULL,
                                           label = "GWAS") {
  declared_genome <- normalize_genome_build(declared_genome)
  values <- NULL
  source <- "not_available"
  if (!is.null(genome_column)) {
    if (!genome_column %in% colnames(table)) {
      stop(label, " genome metadata column not found: ", genome_column, call. = FALSE)
    }
    values <- unique(as.character(table[[genome_column]]))
    source <- paste0("column:", genome_column)
  } else {
    for (attribute in c("genome", "genome_build")) {
      value <- attr(table, attribute, exact = TRUE)
      if (!is.null(value)) {
        values <- unique(as.character(value))
        source <- paste0("attribute:", attribute)
        break
      }
    }
  }
  values <- values[!is.na(values) & nzchar(values)]
  if (!length(values)) {
    return(data.frame(
      input = label, declared_genome = declared_genome,
      stored_genome = NA_character_, metadata_status = "not_available",
      metadata_source = source, stringsAsFactors = FALSE
    ))
  }
  builds <- unique(vapply(values, normalize_genome_build, character(1)))
  if (length(builds) != 1L || builds != declared_genome) {
    stop(
      label, " is declared as ", declared_genome,
      " but its genome metadata reports ", paste(builds, collapse = ", "), ".",
      call. = FALSE
    )
  }
  data.frame(
    input = label, declared_genome = declared_genome,
    stored_genome = builds, metadata_status = "matched",
    metadata_source = source, stringsAsFactors = FALSE
  )
}

chain_config_key <- function(source_genome, target_genome) {
  paste0(normalize_genome_build(source_genome), "_to_", normalize_genome_build(target_genome))
}

resolve_chain_file <- function(config, source_genome, target_genome, must_work = TRUE) {
  source_genome <- normalize_genome_build(source_genome)
  target_genome <- normalize_genome_build(target_genome)
  if (source_genome == target_genome) return(NULL)
  key <- chain_config_key(source_genome, target_genome)
  path <- config$genome$chains[[key]] %||% NULL
  if (is.null(path) || !is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop(
      "Genome conversion ", source_genome, " -> ", target_genome,
      " requires genome.chains.", key, " in the YAML configuration.", call. = FALSE
    )
  }
  path <- resolve_config_path(config, path)
  if (must_work) assert_file(path, paste0("LiftOver chain ", source_genome, " -> ", target_genome)) else path
}

liftover_qc_row <- function(source_genome, target_genome, n_input, n_mapped_once,
                            n_unmapped, n_multimapped, n_duplicate_target, n_retained,
                            chain_file = NULL) {
  data.frame(
    source_genome = source_genome,
    target_genome = target_genome,
    n_input = as.integer(n_input),
    n_mapped_once = as.integer(n_mapped_once),
    n_unmapped = as.integer(n_unmapped),
    n_multimapped = as.integer(n_multimapped),
    n_duplicate_target = as.integer(n_duplicate_target),
    n_retained = as.integer(n_retained),
    chain_file = chain_file %||% NA_character_,
    chain_md5 = if (is.null(chain_file)) NA_character_ else unname(tools::md5sum(chain_file)),
    stringsAsFactors = FALSE
  )
}

import_liftover_chain <- function(path) {
  require_namespaces("rtracklayer", "hg19/hg38 liftOver")
  if (!grepl("\\.gz$", path, ignore.case = TRUE)) return(rtracklayer::import.chain(path))
  unpacked <- tempfile(pattern = "EPIC2-chain-", fileext = ".chain")
  input <- gzfile(path, open = "rb")
  output <- file(unpacked, open = "wb")
  on.exit({
    try(close(input), silent = TRUE)
    try(close(output), silent = TRUE)
    unlink(unpacked)
  }, add = TRUE)
  repeat {
    block <- readBin(input, what = "raw", n = 1024L * 1024L)
    if (!length(block)) break
    writeBin(block, output)
  }
  close(input)
  close(output)
  rtracklayer::import.chain(unpacked)
}

#' Lift genomic ranges between hg19 and hg38
#'
#' Only ranges with exactly one mapping are retained by default. A directional
#' UCSC chain file must be supplied whenever source and target builds differ.
#'
#' @param ranges A named `GRanges` object.
#' @param source_genome Declared source genome.
#' @param target_genome Requested target genome.
#' @param chain_file Directional UCSC chain file.
#' @param unmapped_action Either `drop` or `error`.
#' @param multimapping_action Either `drop` or `error`.
#' @param duplicate_action How duplicate target intervals are handled.
#' @param label Input label used in diagnostics.
#' @return A list containing lifted ranges, retained input indices, and QC counts.
#' @export
lift_genomic_ranges <- function(ranges, source_genome, target_genome, chain_file = NULL,
                                unmapped_action = "drop", multimapping_action = "drop",
                                duplicate_action = "drop", label = "ranges") {
  require_namespaces(c("GenomicRanges", "IRanges", "S4Vectors"), "genome conversion")
  source_genome <- normalize_genome_build(source_genome)
  target_genome <- normalize_genome_build(target_genome)
  unmapped_action <- match.arg(unmapped_action, c("drop", "error"))
  multimapping_action <- match.arg(multimapping_action, c("drop", "error"))
  duplicate_action <- match.arg(duplicate_action, c("drop", "error"))
  if (is.null(names(ranges))) names(ranges) <- paste0("range_", seq_along(ranges))
  if (anyNA(names(ranges)) || any(!nzchar(names(ranges))) || anyDuplicated(names(ranges))) {
    stop(label, " must have unique, non-empty range names before liftOver.", call. = FALSE)
  }
  ranges <- normalize_granges_seqnames(ranges)
  validate_genome_coordinates(
    as.character(GenomicRanges::seqnames(ranges)), GenomicRanges::start(ranges),
    GenomicRanges::end(ranges), source_genome, label
  )

  if (source_genome == target_genome) {
    mapped <- ranges
    input_index <- seq_along(ranges)
    chain_path <- NULL
    n_unmapped <- 0L
    n_multimapped <- 0L
  } else {
    require_namespaces("rtracklayer", "hg19/hg38 liftOver")
    chain_path <- assert_file(chain_file, paste0("LiftOver chain ", source_genome, " -> ", target_genome))
    chain <- import_liftover_chain(chain_path)
    lifted <- rtracklayer::liftOver(ranges, chain)
    counts <- S4Vectors::elementNROWS(lifted)
    n_unmapped <- sum(counts == 0L)
    n_multimapped <- sum(counts > 1L)
    if (n_unmapped && unmapped_action == "error") {
      stop(label, " has ", n_unmapped, " range(s) with no ", target_genome, " mapping.", call. = FALSE)
    }
    if (n_multimapped && multimapping_action == "error") {
      stop(label, " has ", n_multimapped, " range(s) with multiple ", target_genome, " mappings.", call. = FALSE)
    }
    input_index <- which(counts == 1L)
    mapped <- unlist(lifted[input_index], use.names = FALSE)
    names(mapped) <- names(ranges)[input_index]
    mapped <- normalize_granges_seqnames(mapped)
  }

  if (length(mapped)) {
    validate_genome_coordinates(
      as.character(GenomicRanges::seqnames(mapped)), GenomicRanges::start(mapped),
      GenomicRanges::end(mapped), target_genome, paste0(label, " after liftOver")
    )
  }
  target_id <- paste0(
    canonical_chr(as.character(GenomicRanges::seqnames(mapped))), "-",
    GenomicRanges::start(mapped), "-", GenomicRanges::end(mapped)
  )
  duplicate <- duplicated(target_id) | duplicated(target_id, fromLast = TRUE)
  n_duplicate <- sum(duplicate)
  if (n_duplicate && duplicate_action == "error") {
    stop(label, " produces ", n_duplicate, " duplicated target range(s) after liftOver.", call. = FALSE)
  }
  if (n_duplicate) {
    mapped <- mapped[!duplicate]
    input_index <- input_index[!duplicate]
  }
  if (!length(mapped)) {
    stop(label, " has no unique one-to-one ranges after genome harmonization.", call. = FALSE)
  }
  names(mapped) <- names(ranges)[input_index]
  qc <- liftover_qc_row(
    source_genome, target_genome, length(ranges),
    length(ranges) - n_unmapped - n_multimapped, n_unmapped, n_multimapped,
    n_duplicate, length(mapped), chain_path
  )
  list(ranges = mapped, input_index = input_index, qc = qc)
}

harmonize_gwas_coordinates <- function(gwas, columns, source_genome, target_genome,
                                       chain_file = NULL, genome_column = NULL,
                                       unmapped_action = "drop",
                                       multimapping_action = "drop",
                                       duplicate_action = "drop", gwas_name = NULL) {
  label <- if (is.null(gwas_name)) "GWAS" else paste0("GWAS '", gwas_name, "'")
  source_genome <- normalize_genome_build(source_genome)
  target_genome <- normalize_genome_build(target_genome)
  metadata_qc <- check_table_genome_declaration(gwas, source_genome, genome_column, label)
  chromosome <- canonical_chr(gwas[[columns$chr]])
  position <- suppressWarnings(as.integer(gwas[[columns$pos]]))
  validate_genome_coordinates(chromosome, position, position, source_genome, label)
  if (source_genome == target_genome) {
    target_id <- paste0(chromosome, "-", position)
    duplicate <- duplicated(target_id) | duplicated(target_id, fromLast = TRUE)
    n_duplicate <- sum(duplicate)
    if (n_duplicate && duplicate_action == "error") {
      stop(label, " contains ", n_duplicate, " duplicated target coordinate(s).", call. = FALSE)
    }
    input_index <- if (n_duplicate) which(!duplicate) else seq_len(nrow(gwas))
    if (!length(input_index)) stop(label, " has no unique coordinates after genome validation.", call. = FALSE)
    result <- gwas[input_index, , drop = FALSE]
    result$EPIC2_input_chr <- as.character(gwas[[columns$chr]][input_index])
    result$EPIC2_input_pos <- position[input_index]
    result$EPIC2_input_genome <- source_genome
    result$EPIC2_target_genome <- target_genome
    result[[columns$chr]] <- sub("^chr", "", chromosome[input_index], ignore.case = TRUE)
    result[[columns$pos]] <- position[input_index]
    qc <- liftover_qc_row(
      source_genome, target_genome, nrow(gwas), nrow(gwas), 0L, 0L,
      n_duplicate, length(input_index), NULL
    )
    return(list(
      gwas = result, input_index = input_index, liftover_qc = qc,
      metadata_qc = metadata_qc
    ))
  }

  require_namespaces(c("GenomicRanges", "IRanges"), "GWAS genome harmonization")
  ranges <- GenomicRanges::GRanges(
    seqnames = chromosome,
    ranges = IRanges::IRanges(start = position, end = position),
    strand = "*"
  )
  names(ranges) <- paste0("row_", seq_len(nrow(gwas)))
  lifted <- lift_genomic_ranges(
    ranges, source_genome, target_genome, chain_file,
    unmapped_action, multimapping_action, duplicate_action, label
  )
  result <- gwas[lifted$input_index, , drop = FALSE]
  result$EPIC2_input_chr <- as.character(gwas[[columns$chr]][lifted$input_index])
  result$EPIC2_input_pos <- suppressWarnings(as.integer(gwas[[columns$pos]][lifted$input_index]))
  result$EPIC2_input_genome <- source_genome
  result$EPIC2_target_genome <- target_genome
  result[[columns$chr]] <- sub("^chr", "", canonical_chr(
    as.character(GenomicRanges::seqnames(lifted$ranges))
  ), ignore.case = TRUE)
  result[[columns$pos]] <- GenomicRanges::start(lifted$ranges)
  list(gwas = result, input_index = lifted$input_index, liftover_qc = lifted$qc,
       metadata_qc = metadata_qc)
}
