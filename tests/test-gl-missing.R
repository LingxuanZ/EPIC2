library(EPIC2)

# GL -9 is a valid state, including when the pipeline supplies its default
# missing code explicitly. Scalar text sentinel cells remain missing.
gl <- matrix(c("0,-9,-18", "-18,-9,0", "-9,-9,-9", "-9", NA, "."),
             nrow = 1L, dimnames = list("rs1", paste0("sample", 1:6)))
expected_low <- (1e-9 + 2e-18) / (1 + 1e-9 + 1e-18)
for (result in list(normalize_genotype(gl, format = "gl"),
                    normalize_genotype(gl, format = "gl", missing_values = -9))) {
  stopifnot(
    max(abs(result[1, 1:3] - c(expected_low, 2 - expected_low, 1))) < 1e-12,
    all(is.na(result[1, 4:6]))
  )
}
custom <- matrix(c("0,-99,-198", "-99"), 1L,
                 dimnames = list("rs1", c("called", "missing")))
custom_result <- normalize_genotype(custom, format = "gl", missing_values = -99)
stopifnot(custom_result[1, 1] > 0, is.finite(custom_result[1, 1]),
          is.na(custom_result[1, 2]))

# Numeric and character arrays must agree; true NA states make the genotype
# missing. An equal three-state likelihood vector is called, not a sentinel.
gl_array <- array(NA_real_, c(1L, 4L, 3L),
                  dimnames = list("rs1", paste0("sample", 1:4), c("0", "1", "2")))
gl_array[1, 1, ] <- c(0, -9, -18)
gl_array[1, 2, ] <- c(-18, -9, 0)
gl_array[1, 3, ] <- rep(-9, 3)
gl_array[1, 4, ] <- c(0, NA, -18)
character_array <- array(as.character(gl_array), dim(gl_array), dimnames(gl_array))
state_list <- lapply(seq_len(3L), function(i) {
  matrix(gl_array[, , i], 1L, dimnames = dimnames(gl_array)[1:2])
})
for (input in list(gl_array, character_array, state_list)) {
  result <- normalize_genotype(input, format = "gl")
  stopifnot(max(abs(result[1, 1:3] - c(expected_low, 2 - expected_low, 1))) < 1e-12,
            is.na(result[1, 4]))
}
character_array[1, 4, 2] <- "missing"
character_result <- normalize_genotype(
  character_array, format = "gl", missing_values = c(-9, "missing")
)
stopifnot(is.na(character_result[1, 4]), all(is.finite(character_result[1, 1:3])))

# The dosage/GP/PL sentinel rules remain unchanged.
dosage <- matrix(c(0, -9, 2), 1L, dimnames = list("rs1", paste0("s", 1:3)))
stopifnot(is.na(normalize_genotype(dosage, format = "dosage")[1, 2]))
for (format in c("gp", "pl")) {
  states <- array(c(1, -9, 0), c(1L, 1L, 3L),
                  dimnames = list("rs1", "sample1", c("0", "1", "2")))
  stopifnot(is.na(normalize_genotype(states, format = format)[1, 1]))
}

# Test inherited YAML defaults through the public GWAS stage, not only the
# standalone normalizer. No missing-values override should be needed for GL.
local({
  stage_dir <- tempfile("EPIC2-gl-stage-")
  dir.create(stage_dir)
  on.exit(unlink(stage_dir, recursive = TRUE), add = TRUE)
  panel <- rbind(
    rs1 = c("0,-9,-18", "-9,0,-9", "-18,-9,0"),
    rs2 = c("-18,-9,0", "-9,0,-9", "0,-9,-18")
  )
  colnames(panel) <- paste0("sample", 1:3)
  gwas <- data.frame(rsid = c("rs1", "rs2"), chr = 1L, pos = c(100L, 200L),
                     A1 = "A", A2 = "G", MAF = 0.25, Beta = c(0.2, 0.3),
                     SE = 0.1, P = 0.1)
  saveRDS(panel, file.path(stage_dir, "panel.rds"))
  saveRDS(gwas, file.path(stage_dir, "gwas.rds"))
  yaml::write_yaml(list(
    genome = list(build = "hg19"),
    project = list(output_dir = "results"),
    atac = list(input = "unused.rds"),
    gwas = list(genotype_format = "gl", ld_prune = FALSE,
                traits = list(T = list(summary = "gwas.rds", genotype = "panel.rds")))
  ), file.path(stage_dir, "config.yml"))
  config <- read_EPIC2_config(file.path(stage_dir, "config.yml"))
  stopifnot(config$gwas$genotype_missing_values == -9)
  paths <- run_gwas_stage(config)
  retained <- readRDS(paths$T$gwas)
  retained_panel <- readRDS(paths$T$genotype)
  stopifnot(identical(retained$rsid, c("rs1", "rs2")),
            identical(rownames(retained_panel), retained$rsid),
            all(is.finite(retained_panel)),
            all(retained$EPIC2_reference_missing_rate == 0))
})

message("GL likelihood-state missingness and pipeline-default regressions passed.")
