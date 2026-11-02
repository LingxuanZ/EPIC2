library(EPIC2)

expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NA_character_ }, error = conditionMessage)
  stopifnot(!is.na(message), grepl(pattern, message, fixed = TRUE))
}

# V4: quantiles exclude infinite raw scores, then cap those scores.
input <- data.frame(pval = c(seq(0.02, 0.98, length.out = 95), rep(0, 5), 0, NA),
                    auc = c(rep(0.7, 100), 0.3, 0.7))
actual <- make_one_sided_wrs(input)
p <- ifelse(input$auc >= 0.5, input$pval / 2, 1 - input$pval / 2)
raw <- qnorm(p, lower.tail = FALSE)
limits <- quantile(raw[is.finite(raw)], c(0.02, 0.98), names = FALSE)
expected <- pmin(pmax(raw, limits[1]), limits[2])
stopifnot(isTRUE(all.equal(actual$z_one_sided_enriched_trunc, expected)),
          all(actual$z_trunc_upper == limits[2]),
          all(is.infinite(actual$z_one_sided_enriched[96:101])))
expect_error(make_one_sided_wrs(data.frame(pval = c(0, 0), auc = c(0.7, 0.3))), "No finite")
expect_error(make_one_sided_wrs(data.frame(pval = -0.1, auc = 0.7)), "[0, 1]")

# V4: count threshold uses only cells already passing TSS and FRIP.
metadata <- data.frame(nCount_peaks = c(1, 1, 100, 101, 1000),
                       TSS.enrichment = c(0, 5, 5, 5, 5), FRIP = c(0.8, 0.1, 0.8, 0.8, 0.8))
cutoff <- EPIC2:::atac_count_cutoff(metadata, 0.01, assay = "peaks")
stopifnot(isTRUE(all.equal(cutoff, unname(quantile(c(100, 101, 1000), 0.01)))),
          cutoff > quantile(metadata$nCount_peaks, 0.01))

# V4: all AVAILABLE valid nominal component P-values, equal weights.
component <- function(p) data.frame(trait = "T", group = "G", p = p)
out <- combine_component_tests(component(0.01), component(0.02), component(0.03), component(NA))
expected <- pcauchy(mean(tan((0.5 - c(0.01, 0.02, 0.03)) * pi)), lower.tail = FALSE)
stopifnot(out$n_tests_combined == 3, abs(out$p_cauchy - expected) < 1e-14)
out <- combine_component_tests(component(0.01), component(0.02), component(0.03), component(-1), require_all = TRUE)
stopifnot(out$n_tests_combined == 3, is.na(out$p_cauchy))

# Label parsing must preserve the identity of the counted allele, so signed LD
# is consistent with the GWAS Z-scores in the peak quadratic form.
gwas <- data.frame(rsid = c("rs1", "rs2"), chr = 1L, pos = c(100, 110),
                    A1 = c("A", "T"), A2 = c("G", "C"), MAF = 0.3,
                    Beta = c(0.2, -0.3), SE = 0.1)
labels <- rbind(c("AA", "AG", "GG", "AG"), c("CC", "CT", "TT", "CC"))
dimnames(labels) <- list(gwas$rsid, paste0("s", 1:4))
dosage <- rbind(c(2, 1, 0, 1), c(0, 1, 2, 0))
dimnames(dosage) <- dimnames(labels)
parsed <- harmonize_gwas(gwas, labels, genotype_format = "allele", ld_prune = FALSE, genome_build = "hg19")
direct <- harmonize_gwas(gwas, dosage, genotype_format = "dosage", ld_prune = FALSE, genome_build = "hg19")
stopifnot(isTRUE(all.equal(parsed$genotype, direct$genotype)), parsed$summary$n_allele_dosage_flips == 1)
bad <- labels; bad[1, 1] <- "AC"
expect_error(harmonize_gwas(gwas, bad, genotype_format = "allele", ld_prune = FALSE, genome_build = "hg19"), "more than two")
missing <- labels; missing[1, 1] <- NA_character_
stopifnot(is.na(normalize_genotype(missing, format = "allele")[1, 1]))

# Carrier fraction is not MAF: six heterozygotes among 100 diploid samples.
carrier_case <- matrix(c(rep(1, 6), rep(0, 94)), nrow = 1)
stopifnot(EPIC2:::reference_maf(carrier_case) == 0.03, mean(carrier_case > 0) == 0.06)

# The stated pruning order can differ from removing ambiguous SNPs first.
# A correlated valid neighbor is displaced by the first ambiguous SNP, then
# that ambiguous SNP is removed; the independent SNP survives.
pruning_gwas <- data.frame(rsid = c("ambiguous", "neighbor", "independent"),
                           chr = 1L, pos = c(100, 200, 300), A1 = "A", A2 = c("T", "G", "G"),
                           MAF = 0.3, Beta = 0.2, SE = 0.1)
pruning_dosage <- rbind(rep(c(0, 1, 2), 4), rep(c(0, 1, 2), 4), rep(c(0, 1, 2), each = 4))
dimnames(pruning_dosage) <- list(pruning_gwas$rsid, paste0("s", 1:12))
pruned <- harmonize_gwas(pruning_gwas, pruning_dosage, genotype_format = "dosage",
                         genome_build = "hg19", ld_method = "corr", ld_start_pos = "first")
stopifnot(identical(pruned$gwas$rsid, "independent"), pruned$summary$n_after_ld_prune == 2,
          pruned$summary$n_after_ambiguous == 1)

# Optional S-LDSC comparison tests the positive coefficient, as in V4.
sldsc <- EPIC2:::standardize_sldsc_summary(data.frame(trait = "T", group = c("a", "b"), coef_z = c(2, -2)))
stopifnot(max(abs(sldsc$coefficient_p - pnorm(c(2, -2), lower.tail = FALSE))) < 1e-15)

# A missing observation cannot be handled by slicing an existing whitening
# matrix: the matrix must be rebuilt from the retained covariance.
set.seed(902)
n <- 20
ids <- paste0("f", seq_len(n))
W <- diag(seq(0.7, 1.8, length.out = n)); dimnames(W) <- list(ids, ids)
x <- setNames(rnorm(n), ids); b <- setNames(rnorm(n), ids)
y <- setNames(2 + 0.3*x + 0.2*b + rnorm(n), ids)
fit <- fit_whitened_groups(y, list(G = x), b, W)
X <- W %*% cbind(intercept = 1, x, b)
reference <- summary(lm(as.numeric(W %*% y) ~ X - 1))$coefficients[2, ]
stopifnot(abs(fit$beta - reference[1]) < 1e-12, abs(fit$se - reference[2]) < 1e-12)
x[2] <- NA
bad_fit <- fit_whitened_groups(y, list(G = x), b, W)
stopifnot(is.na(bad_fit$p), grepl("rebuild covariance", bad_fit$status))
message("Manuscript WRS, QC order, ACAT, allele orientation, MAF, and model consistency checks passed.")
