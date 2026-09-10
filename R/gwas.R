is_strand_ambiguous <- function(a1, a2) {
  pair <- paste0(toupper(as.character(a1)), toupper(as.character(a2)))
  pair %in% c("AT", "TA", "CG", "GC")
}

reference_maf <- function(genotype, ploidy = 2) {
  ploidy <- validate_genotype_ploidy(ploidy)
  dosage <- as.matrix(genotype)
  if (!is.numeric(dosage)) stop("Reference genotypes must be numeric allele dosages.", call. = FALSE)
  if (any(dosage < 0 | dosage > ploidy, na.rm = TRUE)) {
    stop("Reference genotypes must be dosages between 0 and ploidy.", call. = FALSE)
  }
  af <- rowMeans(dosage, na.rm = TRUE) / ploidy
  af[!is.finite(af)] <- NA_real_
  pmin(af, 1 - af)
}

snprelate_pruning_options <- function(threshold, window_bp, method = "composite",
                                      maf_min = 0.005, missing_rate_max = 0.01,
                                      remove_monomorphic = TRUE, slide_max_n = NULL,
                                      start_pos = "random.f500", num_threads = 1L) {
  list(
    autosome.only = TRUE,
    remove.monosnp = remove_monomorphic,
    maf = maf_min,
    missing.rate = missing_rate_max,
    method = method,
    slide.max.bp = as.integer(window_bp),
    slide.max.n = if (is.null(slide_max_n)) NA_integer_ else as.integer(slide_max_n),
    ld.threshold = threshold,
    start.pos = start_pos,
    num.thread = as.integer(num_threads),
    verbose = FALSE
  )
}

ld_prune_snprelate <- function(gwas, genotype, columns, threshold, window_bp, seed, ploidy,
                              noninteger = c("error", "round"),
                              method = "composite", maf_min = 0.005,
                              missing_rate_max = 0.01, remove_monomorphic = TRUE,
                              slide_max_n = NULL, start_pos = "random.f500",
                              num_threads = 1L,
                              output_dir = tempdir()) {
  require_namespaces("SNPRelate", "GWAS LD pruning")
  noninteger <- match.arg(noninteger)
  pruning_genotype <- as_snprelate_hard_calls(genotype, ploidy, noninteger)
  gds_path <- tempfile(pattern = "EPIC2-pruning-", tmpdir = output_dir, fileext = ".gds")
  on.exit(unlink(gds_path), add = TRUE)

  rsid <- as.character(gwas[[columns$rsid]])
  allele <- paste0(gwas[[columns$a1]], "/", gwas[[columns$a2]])
  SNPRelate::snpgdsCreateGeno(
    gds_path,
    genmat = pruning_genotype,
    sample.id = colnames(genotype),
    snp.id = rsid,
    snp.chromosome = chr_to_order(gwas[[columns$chr]]),
    snp.position = as.integer(gwas[[columns$pos]]),
    snp.allele = allele,
    snpfirstdim = TRUE
  )
  handle <- SNPRelate::snpgdsOpen(gds_path)
  on.exit(try(SNPRelate::snpgdsClose(handle), silent = TRUE), add = TRUE)
  set.seed(seed)
  options <- snprelate_pruning_options(
    threshold = threshold,
    window_bp = window_bp,
    method = method,
    maf_min = maf_min,
    missing_rate_max = missing_rate_max,
    remove_monomorphic = remove_monomorphic,
    slide_max_n = slide_max_n,
    start_pos = start_pos,
    num_threads = num_threads
  )
  selected <- do.call(SNPRelate::snpgdsLDpruning, c(list(gdsobj = handle), options))
  selected <- unlist(selected, use.names = FALSE)
  match(selected, rsid, nomatch = 0L)
}

#' Harmonize and prune one GWAS with a reference genotype panel
#'
#' @param gwas GWAS summary-statistics data frame.
#' @param genotype Variant-by-sample dosage matrix with rsIDs as row names.
#' @param columns Column-name mapping from the configuration.
#' @param maf_min Minimum GWAS and reference-panel minor allele frequency.
#' @param remove_ambiguous Remove A/T and C/G variants.
#' @param ld_prune Perform SNPRelate LD pruning.
#' @param ld_threshold LD pruning threshold.
#' @param ld_window_bp LD pruning window in base pairs.
#' @param seed Random seed passed to SNPRelate.
#' @param genotype_ploidy Maximum genotype dosage.
#' @param genotype_format Genotype encoding accepted by [normalize_genotype()].
#' @param genotype_orientation Location of variants in the input matrix.
#' @param genotype_missing_values Input values representing missing genotypes.
#' @param genotype_variant_id_column Optional variant ID column for data-frame genotypes.
#' @param ld_prune_noninteger How SNPRelate pruning handles non-integer dosage.
#' @param gwas_name Optional GWAS name used in diagnostics and the QC summary.
#' @param source_genome Declared GWAS coordinate build.
#' @param target_genome Requested output coordinate build.
#' @param chain_file Directional UCSC chain file when builds differ.
#' @param genome_column Optional GWAS column containing build metadata to verify.
#' @param liftover_unmapped_action Drop or reject variants without a target mapping.
#' @param liftover_multimapping_action Drop or reject variants with multiple mappings.
#' @param liftover_duplicate_action Drop or reject duplicated target coordinates.
#' @param reference_missing_rate_max Maximum missing-genotype rate per reference variant.
#' @param ld_prune_tool LD-pruning implementation; currently `"snprelate"`.
#' @param ld_method SNPRelate LD statistic: `"composite"`, `"r"`, `"dprime"`, or `"corr"`.
#' @param ld_slide_max_n Optional maximum number of variants in a pruning window.
#' @param ld_prune_maf_min SNPRelate's internal MAF threshold during pruning.
#' @param ld_remove_monomorphic Ask SNPRelate to remove monomorphic variants.
#' @param ld_start_pos SNPRelate window-start strategy.
#' @param ld_num_threads Number of SNPRelate pruning threads.
#' @return A list containing aligned `gwas`, `genotype`, and a QC summary.
#' @export
harmonize_gwas <- function(gwas, genotype, columns = EPIC2_defaults()$gwas$columns,
                           maf_min = 0.05, remove_ambiguous = TRUE, ld_prune = TRUE,
                           ld_threshold = 0.8, ld_window_bp = 500000L, seed = 1000L,
                           genotype_ploidy = 2, genotype_format = "auto",
                           genotype_orientation = "auto", genotype_missing_values = -9,
                           genotype_variant_id_column = NULL,
                           ld_prune_noninteger = "error", gwas_name = NULL,
                           source_genome = "hg19", target_genome = source_genome,
                           chain_file = NULL, genome_column = NULL,
                           liftover_unmapped_action = "drop",
                           liftover_multimapping_action = "drop",
                           liftover_duplicate_action = "drop",
                           reference_missing_rate_max = 0.01,
                           ld_prune_tool = "snprelate", ld_method = "composite",
                           ld_slide_max_n = NULL, ld_prune_maf_min = 0.005,
                           ld_remove_monomorphic = TRUE, ld_start_pos = "random.f500",
                           ld_num_threads = 1L) {
  setting_errors <- gwas_setting_errors(list(
    genotype_format = genotype_format,
    genotype_orientation = genotype_orientation,
    genotype_ploidy = genotype_ploidy,
    genotype_variant_id_column = genotype_variant_id_column,
    genome = source_genome,
    genome_column = genome_column,
    maf_min = maf_min,
    reference_missing_rate_max = reference_missing_rate_max,
    remove_ambiguous = remove_ambiguous,
    ld_prune = ld_prune,
    ld_prune_tool = ld_prune_tool,
    ld_prune_noninteger = ld_prune_noninteger,
    ld_method = ld_method,
    ld_threshold = ld_threshold,
    ld_window_bp = ld_window_bp,
    ld_slide_max_n = ld_slide_max_n,
    ld_prune_maf_min = ld_prune_maf_min,
    ld_remove_monomorphic = ld_remove_monomorphic,
    ld_start_pos = ld_start_pos,
    ld_num_threads = ld_num_threads
  ), "harmonize_gwas")
  if (length(setting_errors)) {
    stop("Invalid GWAS settings:\n- ", paste(setting_errors, collapse = "\n- "), call. = FALSE)
  }
  source_genome <- normalize_genome_build(source_genome)
  target_genome <- normalize_genome_build(target_genome)
  if (!is.data.frame(gwas)) gwas <- as.data.frame(gwas)
  required <- unlist(columns[c("rsid", "chr", "pos", "a1", "a2", "maf")], use.names = FALSE)
  missing <- setdiff(required, colnames(gwas))
  if (length(missing)) stop("Missing GWAS columns: ", paste(missing, collapse = ", "), call. = FALSE)
  input_label <- if (is.null(gwas_name)) "GWAS" else paste0("GWAS '", gwas_name, "'")
  input_genome_metadata_qc <- check_table_genome_declaration(
    gwas, source_genome, genome_column, input_label
  )
  validate_genome_coordinates(
    gwas[[columns$chr]], gwas[[columns$pos]], gwas[[columns$pos]],
    source_genome, input_label
  )
  genotype_ploidy <- validate_genotype_ploidy(genotype_ploidy)
  genotype <- normalize_genotype(
    genotype,
    variant_ids = as.character(gwas[[columns$rsid]]),
    format = genotype_format,
    ploidy = genotype_ploidy,
    orientation = genotype_orientation,
    missing_values = genotype_missing_values,
    variant_id_column = genotype_variant_id_column
  )
  normalized_format <- attr(genotype, "genotype_format")
  validate_named_matrix(genotype, "reference genotype")

  input_n <- nrow(gwas)
  rsid <- as.character(gwas[[columns$rsid]])
  keep <- !is.na(rsid) & nzchar(rsid) & !duplicated(rsid) & rsid %in% rownames(genotype)
  gwas <- gwas[keep, , drop = FALSE]
  genotype <- genotype[match(as.character(gwas[[columns$rsid]]), rownames(genotype)), , drop = FALSE]
  stopifnot(identical(as.character(gwas[[columns$rsid]]), rownames(genotype)))
  after_reference_match <- nrow(gwas)

  reference_missing_rate <- rowMeans(is.na(genotype))
  keep <- is.finite(reference_missing_rate) &
    reference_missing_rate <= reference_missing_rate_max
  gwas <- gwas[keep, , drop = FALSE]
  genotype <- genotype[keep, , drop = FALSE]
  reference_missing_rate <- reference_missing_rate[keep]
  gwas$EPIC2_reference_missing_rate <- reference_missing_rate
  after_reference_missing <- nrow(gwas)

  gwas_maf <- suppressWarnings(as.numeric(gwas[[columns$maf]]))
  invalid_gwas_maf <- is.finite(gwas_maf) & (gwas_maf < 0 | gwas_maf > 0.5)
  if (any(invalid_gwas_maf)) {
    label <- if (is.null(gwas_name)) "GWAS" else paste0("GWAS '", gwas_name, "'")
    stop(
      label, " has ", sum(invalid_gwas_maf), " MAF value(s) outside [0, 0.5]. ",
      "If this column is allele frequency rather than MAF, convert it with pmin(AF, 1 - AF).",
      call. = FALSE
    )
  }
  ref_maf <- reference_maf(genotype, genotype_ploidy)
  gwas$EPIC2_reference_maf <- ref_maf
  keep <- is.finite(gwas_maf) & is.finite(ref_maf) & gwas_maf >= maf_min & ref_maf >= maf_min
  gwas <- gwas[keep, , drop = FALSE]
  genotype <- genotype[keep, , drop = FALSE]
  after_maf <- nrow(gwas)

  if (remove_ambiguous) {
    ambiguous <- is_strand_ambiguous(gwas[[columns$a1]], gwas[[columns$a2]])
    ambiguous[is.na(ambiguous)] <- TRUE
    gwas <- gwas[!ambiguous, , drop = FALSE]
    genotype <- genotype[!ambiguous, , drop = FALSE]
  }
  after_ambiguous <- nrow(gwas)

  coordinate_result <- harmonize_gwas_coordinates(
    gwas = gwas,
    columns = columns,
    source_genome = source_genome,
    target_genome = target_genome,
    chain_file = chain_file,
    genome_column = genome_column,
    unmapped_action = liftover_unmapped_action,
    multimapping_action = liftover_multimapping_action,
    duplicate_action = liftover_duplicate_action,
    gwas_name = gwas_name
  )
  gwas <- coordinate_result$gwas
  genotype <- genotype[coordinate_result$input_index, , drop = FALSE]
  after_liftover <- nrow(gwas)

  chr_ord <- chr_to_order(gwas[[columns$chr]])
  pos <- suppressWarnings(as.integer(gwas[[columns$pos]]))
  keep <- chr_ord %in% 1:22 & is.finite(pos)
  gwas <- gwas[keep, , drop = FALSE]
  genotype <- genotype[keep, , drop = FALSE]
  ord <- order(chr_to_order(gwas[[columns$chr]]), as.integer(gwas[[columns$pos]]))
  gwas <- gwas[ord, , drop = FALSE]
  genotype <- genotype[ord, , drop = FALSE]
  after_autosomal <- nrow(gwas)

  if (ld_prune && nrow(gwas) > 1L) {
    selected <- switch(
      ld_prune_tool,
      snprelate = ld_prune_snprelate(
        gwas, genotype, columns, ld_threshold, ld_window_bp, seed,
        genotype_ploidy, ld_prune_noninteger,
        method = ld_method,
        maf_min = ld_prune_maf_min,
        missing_rate_max = reference_missing_rate_max,
        remove_monomorphic = ld_remove_monomorphic,
        slide_max_n = ld_slide_max_n,
        start_pos = ld_start_pos,
        num_threads = ld_num_threads
      ),
      stop("Unsupported LD-pruning tool: ", ld_prune_tool, call. = FALSE)
    )
    selected <- selected[selected > 0L]
    gwas <- gwas[selected, , drop = FALSE]
    genotype <- genotype[selected, , drop = FALSE]
  }
  after_ld_prune <- nrow(gwas)

  rownames(genotype) <- as.character(gwas[[columns$rsid]])
  gwas$coord_id <- coordinate_id(gwas[[columns$chr]], gwas[[columns$pos]])
  if (anyDuplicated(gwas$coord_id)) {
    stop("GWAS contains duplicated genomic coordinates after harmonization.", call. = FALSE)
  }
  beta_se <- get_beta_se_columns(gwas, columns)
  beta <- suppressWarnings(as.numeric(gwas[[beta_se$beta]]))
  se <- suppressWarnings(as.numeric(gwas[[beta_se$se]]))
  gwas$Z <- beta / se
  keep <- is.finite(gwas$Z) & se > 0
  gwas <- gwas[keep, , drop = FALSE]
  genotype <- genotype[keep, , drop = FALSE]
  stopifnot(identical(as.character(gwas[[columns$rsid]]), rownames(genotype)))

  summary <- data.frame(
    gwas_name = gwas_name %||% NA_character_,
    genotype_format = normalized_format,
    genotype_ploidy = genotype_ploidy,
    source_genome = source_genome,
    target_genome = target_genome,
    genome_metadata_status = input_genome_metadata_qc$metadata_status,
    liftover_chain = coordinate_result$liftover_qc$chain_file,
    liftover_chain_md5 = coordinate_result$liftover_qc$chain_md5,
    n_liftover_unmapped = coordinate_result$liftover_qc$n_unmapped,
    n_liftover_multimapped = coordinate_result$liftover_qc$n_multimapped,
    n_liftover_duplicate_target = coordinate_result$liftover_qc$n_duplicate_target,
    n_reference_samples = ncol(genotype),
    maf_min = maf_min,
    reference_missing_rate_max = reference_missing_rate_max,
    remove_ambiguous = remove_ambiguous,
    ld_prune = ld_prune,
    ld_prune_tool = ld_prune_tool,
    ld_prune_noninteger = ld_prune_noninteger,
    ld_method = ld_method,
    ld_threshold = ld_threshold,
    ld_window_bp = as.integer(ld_window_bp),
    ld_slide_max_n = if (is.null(ld_slide_max_n)) NA_integer_ else as.integer(ld_slide_max_n),
    ld_prune_maf_min = ld_prune_maf_min,
    ld_remove_monomorphic = ld_remove_monomorphic,
    ld_start_pos = ld_start_pos,
    ld_num_threads = as.integer(ld_num_threads),
    n_input = input_n,
    n_reference_matched = after_reference_match,
    n_after_reference_missing = after_reference_missing,
    n_after_maf = after_maf,
    n_after_ambiguous = after_ambiguous,
    n_after_liftover = after_liftover,
    n_after_autosomal = after_autosomal,
    n_after_ld_prune = after_ld_prune,
    n_final = nrow(gwas),
    stringsAsFactors = FALSE
  )
  list(gwas = gwas, genotype = genotype, summary = summary)
}

#' GWAS harmonization and LD pruning
#'
#' @param config Parsed EPIC2 configuration.
#' @param traits Optional GWAS-name subset retained for backward compatibility.
#' @param gwas_name Optional name or names under `gwas.traits` to process.
#' @return Named list of output paths.
#' @export
run_gwas_stage <- function(config, traits = NULL, gwas_name = NULL) {
  traits <- resolve_gwas_names(config, traits = traits, gwas_name = gwas_name)
  output_dir <- stage_output_dir(config, "02_gwas")
  summaries <- list()
  outputs <- setNames(vector("list", length(traits)), traits)

  for (trait in traits) {
    spec <- trait_config(config, trait)
    setting <- function(name) spec[[name]] %||% config$gwas[[name]]
    gwas <- load_analysis_input(
      resolve_config_path(config, spec$summary), spec$summary_object %||% NULL
    )
    genotype <- load_analysis_input(
      resolve_config_path(config, spec$genotype), spec$genotype_object %||% NULL
    )
    source_genome <- normalize_genome_build(setting("genome"))
    target_genome <- normalize_genome_build(config$genome$target)
    chain_file <- resolve_chain_file(config, source_genome, target_genome)
    result <- harmonize_gwas(
      gwas = gwas,
      genotype = genotype,
      columns = deep_merge(config$gwas$columns, spec$columns %||% list()),
      maf_min = setting("maf_min"),
      remove_ambiguous = setting("remove_ambiguous"),
      ld_prune = setting("ld_prune"),
      ld_threshold = setting("ld_threshold"),
      ld_window_bp = setting("ld_window_bp"),
      seed = config$parameters$seed,
      genotype_ploidy = setting("genotype_ploidy"),
      genotype_format = setting("genotype_format"),
      genotype_orientation = setting("genotype_orientation"),
      genotype_missing_values = setting("genotype_missing_values"),
      genotype_variant_id_column = setting("genotype_variant_id_column"),
      ld_prune_noninteger = setting("ld_prune_noninteger"),
      gwas_name = trait,
      source_genome = source_genome,
      target_genome = target_genome,
      chain_file = chain_file,
      genome_column = setting("genome_column"),
      liftover_unmapped_action = config$genome$unmapped_action,
      liftover_multimapping_action = config$genome$multimapping_action,
      liftover_duplicate_action = config$genome$duplicate_action,
      reference_missing_rate_max = setting("reference_missing_rate_max"),
      ld_prune_tool = setting("ld_prune_tool"),
      ld_method = setting("ld_method"),
      ld_slide_max_n = setting("ld_slide_max_n"),
      ld_prune_maf_min = setting("ld_prune_maf_min"),
      ld_remove_monomorphic = setting("ld_remove_monomorphic"),
      ld_start_pos = setting("ld_start_pos"),
      ld_num_threads = setting("ld_num_threads")
    )
    gwas_path <- file.path(output_dir, paste0(trait, ".gwas.rds"))
    genotype_path <- file.path(output_dir, paste0(trait, ".genotype.rds"))
    save_rds_atomic(result$gwas, gwas_path)
    save_rds_atomic(result$genotype, genotype_path)
    summaries[[trait]] <- result$summary
    write_table(result$summary, file.path(output_dir, paste0(trait, ".gwas_qc.tsv")))
    outputs[[trait]] <- list(gwas = gwas_path, genotype = genotype_path)
    message("GWAS stage complete for ", trait, ": ", nrow(result$gwas), " variants")
  }
  summary_path <- file.path(output_dir, "gwas_qc_summary.tsv")
  summary_table <- do.call(rbind, summaries)
  if (file.exists(summary_path)) {
    previous <- utils::read.delim(summary_path, stringsAsFactors = FALSE, check.names = FALSE)
    if ("gwas_name" %in% colnames(previous)) {
      previous <- previous[!previous$gwas_name %in% traits, , drop = FALSE]
      shared <- intersect(colnames(previous), colnames(summary_table))
      if (length(shared) == ncol(summary_table) && ncol(previous) == ncol(summary_table)) {
        previous <- previous[, colnames(summary_table), drop = FALSE]
        summary_table <- rbind(previous, summary_table)
      }
    }
  }
  configured_order <- names(config$gwas$traits)
  summary_table <- summary_table[order(match(summary_table$gwas_name, configured_order)), , drop = FALSE]
  rownames(summary_table) <- NULL
  write_table(summary_table, summary_path)
  outputs
}
