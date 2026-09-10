library(EPIC2)

expect_error_message <- function(expression, pattern) {
  message <- tryCatch({
    force(expression)
    NA_character_
  }, error = conditionMessage)
  stopifnot(!is.na(message), grepl(pattern, message, ignore.case = TRUE))
}

# Common genotype encodings normalize to the same variant-by-sample dosage matrix.
expected_dosage <- rbind(
  rs1 = c(0, 1, 2),
  rs2 = c(2, 1, 0)
)
colnames(expected_dosage) <- c("s1", "s2", "s3")
stopifnot(all(normalize_genotype(expected_dosage, variant_ids = rownames(expected_dosage)) == expected_dosage))
stopifnot(all(normalize_genotype(t(expected_dosage), variant_ids = rownames(expected_dosage)) == expected_dosage))

gt_genotype <- rbind(
  rs1 = c("0/0", "0/1", "1|1"),
  rs2 = c("1/1:20", "0|1:15", "0/0:12")
)
stopifnot(all(normalize_genotype(gt_genotype, format = "gt") == expected_dosage))
gt_with_ad <- gt_genotype
gt_with_ad[1, 2] <- "0/1:10,8"
stopifnot(all(normalize_genotype(gt_with_ad) == expected_dosage))

allele_genotype <- rbind(
  rs1 = c("AA", "AB", "BB"),
  rs2 = c("G/G", "A/G", "A/A")
)
stopifnot(all(normalize_genotype(allele_genotype, format = "allele") == expected_dosage))

gp_genotype <- rbind(
  rs1 = c("1,0,0", "0,1,0", "0,0,1"),
  rs2 = c("0,0,1", "0,1,0", "1,0,0")
)
stopifnot(max(abs(normalize_genotype(gp_genotype, format = "gp") - expected_dosage)) < 1e-12)
stopifnot(max(abs(normalize_genotype(gp_genotype) - expected_dosage)) < 1e-12)

pl_genotype <- rbind(
  rs1 = c("0,200,400", "200,0,200", "400,200,0"),
  rs2 = c("400,200,0", "200,0,200", "0,200,400")
)
stopifnot(max(abs(normalize_genotype(pl_genotype, format = "pl") - expected_dosage)) < 1e-12)
expect_error_message(normalize_genotype(pl_genotype), "ambiguous")

# Character probability arrays preserve explicitly missing genotypes.
gp_array <- array(
  as.character(c(1, NA, 0, NA, 0, NA)),
  dim = c(1, 2, 3),
  dimnames = list("rs1", c("s1", "s2"), c("0", "1", "2"))
)
gp_array_dosage <- normalize_genotype(gp_array, format = "gp")
stopifnot(gp_array_dosage["rs1", "s1"] == 0, is.na(gp_array_dosage["rs1", "s2"]))

logical_genotype <- expected_dosage > 0
expect_error_message(
  normalize_genotype(logical_genotype, variant_ids = rownames(logical_genotype)),
  "ambiguous"
)
binary_numeric <- expected_dosage > 0
storage.mode(binary_numeric) <- "numeric"
expect_error_message(normalize_genotype(binary_numeric), "ambiguous")
stopifnot(all(normalize_genotype(binary_numeric, format = "dosage") == binary_numeric))
expect_error_message(
  normalize_genotype(logical_genotype, format = "carrier"),
  "do not contain enough information"
)
stopifnot(all(EPIC2:::as_snprelate_hard_calls(expected_dosage, 2) == expected_dosage))
haploid <- expected_dosage
haploid[haploid == 2] <- 1
stopifnot(all(EPIC2:::as_snprelate_hard_calls(haploid, 1) == 2 * haploid))
continuous <- expected_dosage
continuous[continuous == 1] <- 0.8
expect_error_message(EPIC2:::as_snprelate_hard_calls(continuous, 2), "requires diploid")
stopifnot(all(EPIC2:::as_snprelate_hard_calls(continuous, 2, "round") == expected_dosage))

# Genome names are explicit and coordinate bounds can catch some bad declarations.
stopifnot(
  identical(normalize_genome_build("GRCh37"), "hg19"),
  identical(normalize_genome_build("hg38"), "hg38")
)
validate_genome_coordinates("chr1", 249000000L, genome = "hg19")
expect_error_message(
  validate_genome_coordinates("chr1", 249000000L, genome = "hg38"),
  "incorrectly declared genome build"
)
expect_error_message(normalize_genome_build("GRCh19"), "Unsupported genome build")

# User genotype matrices can be read from delimited files using an ID column.
genotype_tsv <- tempfile(fileext = ".tsv")
utils::write.table(
  data.frame(variant_id = rownames(expected_dosage), expected_dosage, check.names = FALSE),
  genotype_tsv, sep = "\t", quote = FALSE, row.names = FALSE
)
loaded_genotype <- EPIC2:::load_analysis_input(genotype_tsv, label = "test genotype")
loaded_dosage <- normalize_genotype(
  loaded_genotype, variant_ids = rownames(expected_dosage), format = "dosage",
  variant_id_column = "variant_id"
)
stopifnot(all(loaded_dosage == expected_dosage))

# Every SNPRelate pruning option is explicit rather than inherited silently.
pruning_options <- EPIC2:::snprelate_pruning_options(
  threshold = 0.7,
  window_bp = 250000,
  method = "corr",
  maf_min = 0.02,
  missing_rate_max = 0.05,
  remove_monomorphic = FALSE,
  slide_max_n = NULL,
  start_pos = "first",
  num_threads = 3
)
stopifnot(
  identical(pruning_options$method, "corr"),
  pruning_options$ld.threshold == 0.7,
  pruning_options$slide.max.bp == 250000L,
  is.na(pruning_options$slide.max.n),
  pruning_options$maf == 0.02,
  pruning_options$missing.rate == 0.05,
  identical(pruning_options$remove.monosnp, FALSE),
  identical(pruning_options$start.pos, "first"),
  pruning_options$num.thread == 3L
)

# GWAS filtering uses dosage-based MAF and removes ambiguous variants.
gwas <- data.frame(
  rsid = c("rs1", "rs2", "rs3", "rs4"),
  chr = c(1, 1, 1, 2),
  pos = c(100, 200, 300, 100),
  A1 = c("A", "A", "C", "C"),
  A2 = c("G", "T", "A", "T"),
  MAF = c(0.25, 0.25, 0.01, 0.25),
  Beta = c(0.2, 0.1, 0.1, -0.3),
  SE = c(0.1, 0.1, 0.1, 0.1),
  P = c(0.05, 0.2, 0.3, 0.01)
)
genotype <- rbind(
  rs1 = c(0, 0, 1, 2),
  rs2 = c(0, 1, 1, 0),
  rs3 = c(0, 0, 0, 1),
  rs4 = c(0, 1, 1, 0)
)
colnames(genotype) <- paste0("sample", 1:4)
harmonized <- harmonize_gwas(gwas, genotype, ld_prune = FALSE)
stopifnot(identical(harmonized$gwas$rsid, c("rs1", "rs4")))
stopifnot(identical(rownames(harmonized$genotype), c("rs1", "rs4")))
stopifnot(all.equal(harmonized$gwas$EPIC2_reference_maf, c(0.375, 0.25)))
stopifnot(all(harmonized$gwas$EPIC2_reference_missing_rate == 0))
stopifnot(
  identical(harmonized$summary$source_genome, "hg19"),
  identical(harmonized$summary$target_genome, "hg19"),
  all(harmonized$gwas$EPIC2_input_pos == harmonized$gwas$pos)
)

gwas_with_build <- gwas
gwas_with_build$build <- "GRCh38"
expect_error_message(
  harmonize_gwas(
    gwas_with_build, genotype, ld_prune = FALSE,
    source_genome = "hg19", genome_column = "build"
  ),
  "metadata reports hg38"
)

genotype_with_missing <- genotype
genotype_with_missing["rs1", "sample1"] <- NA_real_
strict_missing <- harmonize_gwas(
  gwas, genotype_with_missing,
  ld_prune = FALSE,
  reference_missing_rate_max = 0.2
)
relaxed_missing <- harmonize_gwas(
  gwas, genotype_with_missing,
  ld_prune = FALSE,
  reference_missing_rate_max = 0.25
)
stopifnot(identical(strict_missing$gwas$rsid, "rs4"))
stopifnot(identical(relaxed_missing$gwas$rsid, c("rs1", "rs4")))

# Each configured GWAS can run independently by name without erasing prior QC rows.
stage_dir <- tempfile("EPIC2-gwas-stage-")
dir.create(stage_dir)
saveRDS(gwas, file.path(stage_dir, "HDL.rds"))
saveRDS(gwas, file.path(stage_dir, "LDL.rds"))
saveRDS(genotype, file.path(stage_dir, "dosage.rds"))
saveRDS(apply(genotype, c(1, 2), function(value) c("0/0", "0/1", "1/1")[[value + 1L]]),
        file.path(stage_dir, "gt.rds"))
config <- EPIC2:::EPIC2_defaults()
config$atac$input <- "unused.rds"
config$.config_dir <- stage_dir
config$output_dir <- file.path(stage_dir, "results")
config$gwas$ld_prune <- FALSE
config$gwas$traits <- list(
  HDL = list(
    summary = "HDL.rds", genotype = "dosage.rds",
    maf_min = 0.2, reference_missing_rate_max = 0.2, ld_method = "corr"
  ),
  LDL = list(summary = "LDL.rds", genotype = "gt.rds", genotype_format = "gt")
)
validate_EPIC2_config(config)
hdl_output <- run_gwas_stage(config, gwas_name = "HDL")
stopifnot(identical(names(hdl_output), "HDL"))
stopifnot(file.exists(file.path(config$output_dir, "02_gwas", "HDL.gwas_qc.tsv")))
stopifnot(!file.exists(file.path(config$output_dir, "02_gwas", "LDL.gwas.rds")))
ldl_output <- run_gwas_stage(config, gwas_name = "LDL")
stopifnot(identical(names(ldl_output), "LDL"))
stage_summary <- utils::read.delim(
  file.path(config$output_dir, "02_gwas", "gwas_qc_summary.tsv"),
  stringsAsFactors = FALSE
)
stopifnot(identical(stage_summary$gwas_name, c("HDL", "LDL")))
hdl_qc <- stage_summary[stage_summary$gwas_name == "HDL", , drop = FALSE]
stopifnot(
  hdl_qc$maf_min == 0.2,
  hdl_qc$reference_missing_rate_max == 0.2,
  identical(hdl_qc$ld_method, "corr")
)
expect_error_message(run_gwas_stage(config, gwas_name = "UNKNOWN"), "Unknown GWAS")

invalid_config <- config
invalid_config$gwas$traits$HDL$ld_method <- "not-a-method"
expect_error_message(validate_EPIC2_config(invalid_config), "ld_method")
invalid_config <- config
invalid_config$atac$qc$frip_min <- 1.1
expect_error_message(validate_EPIC2_config(invalid_config), "frip_min")

hg38_config <- config
hg38_config$genome$target <- "GRCh38"
hg38_config$atac$genome <- "hg38"
hg38_config$gwas$genome <- "GRCh38"
validate_EPIC2_config(hg38_config)

mixed_config <- config
mixed_config$genome$target <- "hg38"
expect_error_message(validate_EPIC2_config(mixed_config), "genome.chains.hg19_to_hg38")

# Ridge inverse square root whitens matrix + ridge * I.
sigma <- matrix(c(1, 0.3, 0.3, 1), 2, 2, dimnames = list(c("a", "b"), c("a", "b")))
inverse <- compute_inverse_sqrt(sigma, ridge = 0.01)
whitened <- inverse %*% (sigma + diag(0.01, 2)) %*% inverse
stopifnot(max(abs(whitened - diag(2))) < 1e-8)

# SNP covariance is symmetric, unit-diagonal, and zero outside the local window.
ld_ids <- c("chr1-100-101", "chr1-200-201", "chr1-2000-2001", "chr2-100-101")
ld_genotype <- rbind(
  c(0, 0, 1, 1, 2, 2),
  c(0, 1, 0, 1, 2, 2),
  c(2, 1, 2, 1, 0, 0),
  c(0, 1, 2, 0, 1, 2)
)
rownames(ld_genotype) <- ld_ids
colnames(ld_genotype) <- paste0("s", 1:6)
ld_gwas <- data.frame(
  coord_id = ld_ids,
  chr = c(1, 1, 1, 2),
  pos = c(100, 200, 2000, 100)
)
ld_sigma <- compute_snp_ld_covariance(ld_genotype, ld_gwas, window_bp = 500)
stopifnot(max(abs(as.matrix(ld_sigma) - t(as.matrix(ld_sigma)))) < 1e-12)
stopifnot(all(diag(as.matrix(ld_sigma)) == 1))
stopifnot(ld_sigma[1, 3] == 0, ld_sigma[1, 4] == 0)
ld_sigma_single <- compute_snp_ld_covariance(
  ld_genotype, ld_gwas, window_bp = 500, correlation_block_size = 1L
)
ld_sigma_wide <- compute_snp_ld_covariance(
  ld_genotype, ld_gwas, window_bp = 500, correlation_block_size = 128L
)
stopifnot(max(abs(as.matrix(ld_sigma_single) - as.matrix(ld_sigma_wide))) < 1e-12)

# Sparse group aggregation is equivalent to direct scaled row sums.
feature_matrix <- Matrix::Matrix(
  matrix(c(1, 0, 2, 1, 0, 3, 0, 1), nrow = 2,
         dimnames = list(c("f1", "f2"), paste0("c", 1:4))),
  sparse = TRUE
)
feature_groups <- setNames(c("A", "A", "B", "B"), colnames(feature_matrix))
feature_scale <- setNames(c(1, 2, 0.5, 1.5), colnames(feature_matrix))
aggregated <- build_group_accessibility(feature_matrix, feature_groups, feature_scale)
manual_a <- Matrix::rowSums(
  feature_matrix[, feature_groups == "A", drop = FALSE] %*%
    Matrix::Diagonal(x = feature_scale[feature_groups == "A"])
)
manual_b <- Matrix::rowSums(
  feature_matrix[, feature_groups == "B", drop = FALSE] %*%
    Matrix::Diagonal(x = feature_scale[feature_groups == "B"])
)
stopifnot(
  max(abs(aggregated$by_group$A - manual_a)) < 1e-12,
  max(abs(aggregated$by_group$B - manual_b)) < 1e-12,
  max(abs(aggregated$baseline - rowMeans(cbind(manual_a, manual_b)))) < 1e-12
)

# The model reproduces the notebook by discarding the transformed intercept and
# allowing an ordinary lm() intercept with the transformed predictors.
ids <- paste0("v", 1:8)
covariance <- outer(seq_along(ids), seq_along(ids), function(i, j) 0.4^abs(i - j))
dimnames(covariance) <- list(ids, ids)
whitening <- compute_inverse_sqrt(covariance, ridge = 0)
predictor <- setNames(c(-2, -1, 0, 1, 2, 3, 4, 5), ids)
baseline <- setNames(c(1, 0, 1, 0, 1, 0, 1, 0), ids)
response <- setNames(2 + 0.7 * predictor + 0.2 * baseline + c(0.1, -0.1, 0.05, -0.05, 0.08, -0.08, 0.03, -0.03), ids)
fit <- fit_whitened_groups(response, list(group1 = predictor), baseline, whitening, "trait")
design <- cbind(intercept = 1, predictor = predictor, baseline = baseline)
transformed_design <- as.matrix(whitening %*% design)
manual_data <- data.frame(
  response = as.numeric(whitening %*% response),
  predictor = transformed_design[, "predictor"],
  baseline = transformed_design[, "baseline"]
)
manual <- stats::lm(response ~ predictor + baseline, data = manual_data)
manual_coefficient <- summary(manual)$coefficients["predictor", ]
fully_transformed <- stats::lm.fit(transformed_design, manual_data$response)
stopifnot(
  abs(fit$beta - manual_coefficient[["Estimate"]]) < 1e-10,
  abs(fit$se - manual_coefficient[["Std. Error"]]) < 1e-10,
  abs(fit$t - manual_coefficient[["t value"]]) < 1e-10,
  abs(fit$p_two_sided - manual_coefficient[["Pr(>|t|)"]]) < 1e-10,
  abs(fit$beta - fully_transformed$coefficients[["predictor"]]) > 1e-6
)

# ACAT returns one result per complete trait/group combination.
component <- function(p) data.frame(trait = "T", group = "G", p = p)
combined <- combine_component_tests(component(0.01), component(0.02), component(0.03), component(0.04))
stopifnot(combined$n_tests_combined == 4L)
stopifnot(is.finite(combined$p_cauchy), combined$p_cauchy > 0, combined$p_cauchy < 0.05)

# Repeated WRS stability reports every repeat pair and feature frequency.
repeated <- do.call(rbind, lapply(1:3, function(index) {
  x <- data.frame(
    target_group = "G",
    feature = c("a", "b", "c"),
    auc = c(0.60, 0.56, 0.50) + c(0, index * 0.001, 0),
    pval_one_sided_enriched = c(0.01, 0.03, 0.5),
    z_one_sided_enriched = c(2.3, 1.9, 0),
    stringsAsFactors = FALSE
  )
  x[["repeat"]] <- index
  x[, c("target_group", "repeat", "feature", "auc", "pval_one_sided_enriched", "z_one_sided_enriched")]
}))
repeated$cells_per_side <- 500L
stability <- summarize_wrs_stability(
  repeated,
  top_n = 2,
  auc_cutoff = 0.55,
  stability_frequency_cutoff = 0.8
)
stopifnot(nrow(stability$pairwise) == 3L)
stopifnot(nrow(stability$feature) == 3L)
stopifnot(stability$feature$selected_frequency[stability$feature$feature == "a"] == 1)
stopifnot(stability$feature$stable[stability$feature$feature == "a"])
stopifnot(!stability$feature$stable[stability$feature$feature == "c"])
balanced_summary <- EPIC2:::summarize_balanced_wrs(repeated, auc_cutoff = 0.60)
stopifnot(all(balanced_summary$auc_cutoff == 0.60))
stopifnot(
  balanced_summary$proportion_auc_ge_cutoff[balanced_summary$feature == "a"] == 1,
  balanced_summary$proportion_auc_ge_cutoff[balanced_summary$feature == "b"] == 0
)

invalid_plot_config <- config
invalid_plot_config$plots$formats <- "jpg"
expect_error_message(validate_EPIC2_config(invalid_plot_config), "plots.formats")
