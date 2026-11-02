library(EPIC2)

# Validate effect estimates and standard errors themselves, not just Beta / SE.
ids <- c("bad_beta_na", "bad_beta_inf", "bad_se_inf", "bad_se_zero",
         "bad_se_negative", "bad_se_na", "bad_z_overflow", "good_a", "good_b")
gwas <- data.frame(
  rsid = ids, chr = 1L, pos = seq_along(ids) * 100L,
  A1 = "A", A2 = "G", MAF = 0.25,
  Beta = c(NA, Inf, 1, 1, 1, 1, 1e308, 0.2, -0.3),
  SE = c(0.1, 0.1, Inf, 0, -0.1, NA, 1e-308, 0.1, 0.1), P = 0.1
)
genotype <- matrix(rep(c(0, 1, 2), each = length(ids)), nrow = length(ids),
                   dimnames = list(ids, paste0("sample", 1:3)))
qc <- harmonize_gwas(gwas, genotype, ld_prune = FALSE, genome_build = "hg19")
stopifnot(
  identical(qc$gwas$rsid, c("good_a", "good_b")),
  identical(rownames(qc$genotype), qc$gwas$rsid),
  max(abs(qc$gwas$Z - c(2, -3))) < 1e-12,
  qc$summary$n_after_autosomal == 9L,
  qc$summary$n_invalid_effect_statistics == 7L,
  qc$summary$n_after_effect_qc == 2L,
  qc$summary$n_after_ld_prune == 2L,
  qc$summary$n_final == 2L
)

# If every effect statistic fails QC, stop before pruning or coordinate output.
# Override only the copied function's pruning dependency, not the namespace.
local({
  pruning_called <- FALSE
  harmonize <- harmonize_gwas
  environment(harmonize) <- list2env(list(
    ld_prune_snprelate = function(...) {
      pruning_called <<- TRUE
      stop("Unexpected pruning call")
    }
  ), parent = environment(harmonize_gwas))
  failure <- tryCatch({
    harmonize(gwas[1:7, ], genotype[1:7, , drop = FALSE],
              ld_prune = TRUE, genome_build = "hg19")
    NA_character_
  }, error = conditionMessage)
  stopifnot(identical(failure, "No variants remain after effect-statistic QC."),
            !pruning_called)
})

# Exercise the actual pruning engine: an invalid leading SNP must never remove
# its valid correlated neighbour. A third SNP keeps pruning active after QC.
if (!requireNamespace("SNPRelate", quietly = TRUE)) {
  message("Skipping pre-pruning regression: SNPRelate is unavailable.")
} else {
  pruning_ids <- c("invalid_first", "valid_correlated", "valid_independent")
  pruning_gwas <- data.frame(
    rsid = pruning_ids, chr = 1L, pos = c(100L, 200L, 300L),
    A1 = "A", A2 = "G", MAF = 0.25,
    Beta = c(NA, 0.2, 0.3), SE = 0.1, P = 0.1
  )
  pruning_genotype <- rbind(
    rep(c(0, 1, 2), 4), rep(c(0, 1, 2), 4), rep(c(0, 1, 2), each = 4)
  )
  dimnames(pruning_genotype) <- list(pruning_ids, paste0("sample", 1:12))
  prune <- function(table, panel) harmonize_gwas(
    table, panel, genome_build = "hg19", ld_prune = TRUE,
    ld_method = "corr", ld_start_pos = "first", ld_threshold = 0.8
  )
  # Establish that pruning would select the first SNP from the correlated pair.
  valid_control <- pruning_gwas
  valid_control$Beta[1] <- 0.1
  control <- prune(valid_control, pruning_genotype)
  stopifnot(
    "invalid_first" %in% control$gwas$rsid,
    !"valid_correlated" %in% control$gwas$rsid
  )

  actual <- prune(pruning_gwas, pruning_genotype)
  expected <- prune(pruning_gwas[-1, ], pruning_genotype[-1, , drop = FALSE])
  stopifnot(
    identical(actual$gwas$rsid, c("valid_correlated", "valid_independent")),
    identical(actual$gwas$rsid, expected$gwas$rsid),
    identical(actual$genotype, expected$genotype),
    actual$summary$n_invalid_effect_statistics == 1L,
    actual$summary$n_after_effect_qc == 2L,
    actual$summary$n_after_ld_prune == 2L
  )
}

# A one-trait rerun must preserve older QC rows whose schema predates the new
# effect-statistic fields, leaving unavailable historical counts as NA.
local({
  stage_dir <- tempfile("EPIC2-qc-schema-")
  dir.create(stage_dir)
  on.exit(unlink(stage_dir, recursive = TRUE), add = TRUE)
  saveRDS(gwas, file.path(stage_dir, "gwas.rds"))
  saveRDS(genotype, file.path(stage_dir, "genotype.rds"))
  spec <- list(summary = "gwas.rds", genotype = "genotype.rds")
  yaml::write_yaml(list(
    genome = list(build = "hg19"), project = list(output_dir = "results"),
    atac = list(input = "unused.rds"),
    gwas = list(ld_prune = FALSE, traits = list(A = spec, B = spec))
  ), file.path(stage_dir, "config.yml"))
  config <- read_EPIC2_config(file.path(stage_dir, "config.yml"))
  run_gwas_stage(config, gwas_name = "A")
  summary_path <- file.path(config$output_dir, "02_gwas", "gwas_qc_summary.tsv")
  old <- utils::read.delim(summary_path, stringsAsFactors = FALSE)
  old$n_invalid_effect_statistics <- NULL
  old$n_after_effect_qc <- NULL
  old$legacy_note <- "earlier-run"
  utils::write.table(old, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

  run_gwas_stage(config, gwas_name = "B")
  merged <- utils::read.delim(summary_path, stringsAsFactors = FALSE)
  stopifnot(
    identical(merged$gwas_name, c("A", "B")),
    merged$n_final[1] == old$n_final,
    is.na(merged$n_invalid_effect_statistics[1]),
    is.na(merged$n_after_effect_qc[1]),
    merged$n_invalid_effect_statistics[2] == 7L,
    merged$n_after_effect_qc[2] == 2L,
    identical(merged$legacy_note[1], "earlier-run"),
    is.na(merged$legacy_note[2])
  )
  # Rerunning both traits also handles an empty set of older rows.
  run_gwas_stage(config)
  replaced <- utils::read.delim(summary_path, stringsAsFactors = FALSE)
  stopifnot(identical(replaced$gwas_name, c("A", "B")),
            all(replaced$n_invalid_effect_statistics == 7L),
            all(is.na(replaced$legacy_note)))
})

message("GWAS effect-statistic QC and pre-pruning regressions passed.")
