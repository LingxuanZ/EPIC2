library(EPIC2)

expect_id_error <- function(expression, pattern) {
  message <- tryCatch({ force(expression); NA_character_ }, error = conditionMessage)
  stopifnot(!is.na(message), grepl(pattern, message, ignore.case = TRUE))
}

dosage <- rbind(rsA = c(0, 1, 2, 1), rsB = c(2, 0, 1, 1))
colnames(dosage) <- paste0("sample", seq_len(ncol(dosage)))
anonymous <- dosage
rownames(anonymous) <- NULL
for (orientation in c("auto", "variants_by_samples")) {
  expect_id_error(
    normalize_genotype(anonymous, variant_ids = rownames(dosage), format = "dosage",
                       orientation = orientation),
    "Genotype variants must have row names"
  )
}
expect_id_error(
  normalize_genotype(t(anonymous), variant_ids = rownames(dosage), format = "dosage",
                     orientation = "samples_by_variants"),
  "Genotype variants must have row names"
)

# Probability arrays and lists must not acquire unverified GWAS IDs either.
gp <- array(0, dim = c(nrow(dosage), ncol(dosage), 3L),
            dimnames = c(dimnames(dosage), list(c("0", "1", "2"))))
for (state in 0:2) gp[, , state + 1L] <- dosage == state
gp_unnamed <- gp
dimnames(gp_unnamed) <- list(NULL, colnames(dosage), c("0", "1", "2"))
expect_id_error(
  normalize_genotype(gp_unnamed, variant_ids = rownames(dosage), format = "gp",
                     orientation = "variants_by_samples"),
  "Genotype variants must have row names"
)
gp_list <- lapply(seq_len(3L), function(state) gp_unnamed[, , state])
expect_id_error(
  normalize_genotype(gp_list, variant_ids = rownames(dosage), format = "gp",
                     orientation = "variants_by_samples"),
  "Genotype variants must have row names"
)
stopifnot(all(normalize_genotype(gp, format = "gp") == dosage))

# Metadata IDs must be removed before automatic genotype-format detection.
dosage_table <- data.frame(rsid = rownames(dosage), dosage, check.names = FALSE)
dosage_result <- normalize_genotype(dosage_table, variant_ids = rownames(dosage),
                                    variant_id_column = "rsid")
gp_strings <- matrix(c("1,0,0", "0,0,1", "0,1,0", "1,0,0",
                       "0,0,1", "0,1,0", "0,1,0", "0,1,0"),
                     nrow = 2L, dimnames = dimnames(dosage))
gp_table <- data.frame(rsid = rownames(dosage), gp_strings, check.names = FALSE)
gp_result <- normalize_genotype(gp_table, variant_ids = rownames(dosage),
                                variant_id_column = "rsid")
stopifnot(
  identical(dimnames(dosage_result), dimnames(dosage)),
  identical(attr(dosage_result, "genotype_format"), "dosage"),
  identical(attr(gp_result, "genotype_format"), "gp"),
  all(dosage_result == dosage), all(gp_result == dosage)
)

# Reference-only SNPs are allowed; missing reference SNPs are skipped by ID.
# The surviving rows must remain aligned after coordinate sorting.
gwas <- data.frame(rsid = c("rsB", "rsMissing", "rsA"), chr = 1L,
                   pos = c(200L, 300L, 100L), A1 = "A", A2 = "G",
                   MAF = 0.25, Beta = c(0.2, 0.3, 0.4), SE = 0.1, P = 0.01)
reference <- rbind(rsExtra = c(1, 2, 0, 1), dosage[c("rsB", "rsA"), ])
for (input in list(reference, t(reference))) {
  result <- harmonize_gwas(gwas, input, ld_prune = FALSE, genome_build = "hg19")
  stopifnot(
    identical(result$gwas$rsid, c("rsA", "rsB")),
    identical(rownames(result$genotype), result$gwas$rsid),
    all(result$genotype == dosage),
    result$summary$n_input == 3L, result$summary$n_reference_matched == 2L,
    result$summary$n_final == 2L
  )
}
anonymous_gwas <- gwas[c(3L, 1L), ]
expect_id_error(
  harmonize_gwas(anonymous_gwas, anonymous, ld_prune = FALSE, genome_build = "hg19"),
  "Genotype variants must have row names"
)
duplicate <- dosage
rownames(duplicate) <- c("rsA", "rsA")
expect_id_error(
  normalize_genotype(duplicate, format = "dosage"), "non-missing and unique"
)

message("Genotype identifier and ID-column detection regressions passed.")
