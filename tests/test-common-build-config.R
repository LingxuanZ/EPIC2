# Run with Rscript tests/test-common-build-config.R from any working directory,
# or source this file. R CMD check uses the installed namespace when source R/
# is unavailable. No real analysis is required.
test_source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(test_source_file)) {
  script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_argument)) test_source_file <- sub("^--file=", "", script_argument[[1L]])
}
test_package_dir <- Sys.getenv("EPIC2_TEST_SOURCE_DIR", unset = "")
if (!nzchar(test_package_dir)) {
  test_package_dir <- if (is.null(test_source_file)) {
    "."
  } else {
    file.path(dirname(test_source_file), "..")
  }
}
test_package_dir <- normalizePath(test_package_dir, mustWork = TRUE)

run_common_build_config_tests <- function(package_dir) {
  if (!requireNamespace("yaml", quietly = TRUE)) stop("These configuration tests require yaml.")
  env <- new.env(parent = globalenv())
  source_files <- file.path(package_dir, "R", c("utils.R", "genome.R", "config.R"))
  if (all(file.exists(source_files))) {
    for (path in source_files) sys.source(path, envir = env)
  } else {
    # R CMD check runs tests against the installed namespace, without source R/.
    namespace <- asNamespace("EPIC2")
    for (name in ls(namespace, all.names = TRUE)) {
      assign(name, get(name, envir = namespace, inherits = FALSE), envir = env)
    }
  }
  checks <- 0L
  check <- function(value, label) {
    if (!isTRUE(value)) stop("FAIL: ", label, call. = FALSE)
    checks <<- checks + 1L
    invisible(TRUE)
  }
  expect_error <- function(expression, expected) {
    message <- tryCatch({ force(expression); NA_character_ }, error = conditionMessage)
    check(!is.na(message) && all(vapply(expected, function(text) {
      grepl(text, message, fixed = TRUE)
    }, logical(1))), paste0("expected error: ", paste(expected, collapse = "; "),
                           "; got: ", message))
    invisible(message)
  }
  config_for <- function(build = "hg38") {
    config <- env$EPIC2_defaults()
    config$genome$build <- build
    config$atac$input <- "atac.rds"
    config$gwas$traits <- list(
      HDL = list(summary = "HDL.rds", genotype = "reference.rds"),
      LDL = list(summary = "LDL.rds", genotype = "reference.rds", ld_prune = FALSE)
    )
    config
  }
  set_path <- function(config, path, value) {
    if (length(path) == 1L) {
      config[path] <- list(value)
    } else {
      config[[path[[1L]]]] <- set_path(config[[path[[1L]]]], path[-1L], value)
    }
    config
  }

  defaults <- env$EPIC2_defaults()
  check(is.null(defaults$genome$build) && is.null(defaults$genome$target) &&
          is.null(defaults$atac$genome) && is.null(defaults$gwas$genome),
        "no implicit genome build in defaults")
  check(is.null(defaults$genome$chains) && is.null(defaults$genome$unmapped_action) &&
          is.null(defaults$genome$multimapping_action) &&
          identical(defaults$genome$duplicate_action, "drop"),
        "conversion options optional; duplicate action retained")
  expect_error(env$validate_EPIC2_config(defaults), "genome.build")

  for (declared in c("hg19", "GRCh37", "hg38", "GRCh38")) {
    config <- config_for(declared)
    before <- serialize(config, NULL)
    normalized <- env$validate_EPIC2_config(config)
    canonical <- if (declared %in% c("hg19", "GRCh37")) "hg19" else "hg38"
    check(identical(normalized$genome$build, canonical) &&
            identical(normalized$genome$target, canonical) &&
            identical(normalized$atac$genome, canonical) &&
            identical(normalized$gwas$genome, canonical) &&
            all(vapply(normalized$gwas$traits, function(x) identical(x$genome, canonical), logical(1))),
          paste("canonical aliases:", declared))
    check(identical(serialize(config, NULL), before), "validation does not mutate supplied config")
    check(identical(env$validate_EPIC2_config(normalized), normalized), "normalization is idempotent")
    check(identical(normalized$gwas$traits$LDL$ld_prune, FALSE), "unrelated trait setting preserved")
  }

  legacy_paths <- list(
    c("genome", "target"), c("atac", "genome"), c("gwas", "genome"),
    c("gwas", "traits", "HDL", "genome"), c("gwas", "traits", "LDL", "genome")
  )
  for (path in legacy_paths) {
    alias <- paste(path, collapse = ".")
    config <- set_path(config_for(), path, "GRCh38")
    check(identical(env$validate_EPIC2_config(config)$genome$build, "hg38"),
          paste("same-build legacy alias accepted:", alias))
    config <- set_path(config_for(), path, "hg19")
    expect_error(env$validate_EPIC2_config(config), c(alias, "genome.build", "no automatic"))
    for (invalid in list("", NA_character_, TRUE, 38L, c("hg19", "hg38"), list("hg38"), "mm10")) {
      config <- set_path(config_for(), path, invalid)
      expect_error(env$validate_EPIC2_config(config), paste0("Invalid ", alias))
    }
  }

  config <- config_for()
  for (path in legacy_paths) config <- set_path(config, path, "hg38")
  config$genome$build <- NULL
  expect_error(env$validate_EPIC2_config(config), c("genome.build", "Legacy genome fields"))
  for (invalid in list(NULL, "", NA_character_, TRUE, 38L, c("hg19", "hg38"), list("hg38"), "mm10")) {
    config <- set_path(config_for(), c("genome", "build"), invalid)
    expect_error(env$validate_EPIC2_config(config), "genome.build")
  }
  config <- config_for()
  config$genome <- list(building = "hg38")
  expect_error(env$validate_EPIC2_config(config), "genome.build")

  for (empty in list(NULL, list(), list(hg19_to_hg38 = NULL, hg38_to_hg19 = NULL))) {
    config <- set_path(config_for(), c("genome", "chains"), empty)
    check(identical(env$validate_EPIC2_config(config)$genome$build, "hg38"),
          "empty legacy chain declaration accepted")
  }
  for (nonempty in list("chain.gz", list(hg19_to_hg38 = "missing.chain.gz"),
                       list(unexpected_key = "missing.chain"), list(hg19_to_hg38 = NA_character_))) {
    config <- set_path(config_for(), c("genome", "chains"), nonempty)
    expect_error(env$validate_EPIC2_config(config), c("no automatic", "remove genome.chains"))
  }
  for (option in c("unmapped_action", "multimapping_action")) {
    for (value in list(NULL, "drop", "error")) {
      config <- set_path(config_for(), c("genome", option), value)
      check(identical(env$validate_EPIC2_config(config)$genome$build, "hg38"),
            paste("optional legacy policy:", option))
    }
    config <- set_path(config_for(), c("genome", option), "invalid")
    expect_error(env$validate_EPIC2_config(config), paste0("genome.", option))
  }
  config <- config_for()
  config$genome$duplicate_action <- "error"
  check(identical(env$validate_EPIC2_config(config)$genome$duplicate_action, "error"),
        "duplicate error policy retained")
  config$genome$duplicate_action <- NULL
  expect_error(env$validate_EPIC2_config(config), "genome.duplicate_action")

  # Genome normalization must not bypass the existing non-genome checks.
  config <- config_for(); config$atac$qc$frip_min <- 1.1
  expect_error(env$validate_EPIC2_config(config), "atac.qc.frip_min")
  config <- config_for(); config$ld$ridge <- -1
  expect_error(env$validate_EPIC2_config(config), "ld.ridge")
  config <- config_for(); config$gwas$traits$HDL$ld_method <- "invalid"
  expect_error(env$validate_EPIC2_config(config), "gwas.traits.HDL.ld_method")
  config <- config_for(); config$parameters$workers <- 1.5
  expect_error(env$validate_EPIC2_config(config), "parameters.workers")

  temporary <- tempfile("EPIC2-common-build-config-")
  dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  inputs <- c("atac.rds", "HDL.rds", "reference.rds")
  invisible(file.create(file.path(temporary, inputs)))
  user_config <- list(
    genome = list(build = "GRCh38"), atac = list(input = inputs[[1L]]),
    gwas = list(traits = list(HDL = list(summary = inputs[[2L]], genotype = inputs[[3L]])))
  )
  yaml_file <- file.path(temporary, "config.yml")
  yaml::write_yaml(user_config, yaml_file)
  parsed <- env$read_EPIC2_config(yaml_file, check_files = TRUE)
  check(identical(parsed$genome$build, "hg38") && identical(parsed$atac$genome, "hg38") &&
          identical(parsed$gwas$traits$HDL$genome, "hg38"),
        "YAML global-only config uses returned canonical configuration")
  check(identical(parsed$output_dir, file.path(dirname(normalizePath(yaml_file)), "results")),
        "relative output resolution unchanged")
  check(!dir.exists(parsed$output_dir), "validation does not create an output directory")
  direct <- env$deep_merge(env$EPIC2_defaults(), user_config)
  normalized <- env$validate_EPIC2_config(direct)
  check(identical(normalized$genome, parsed$genome) && identical(normalized$atac, parsed$atac) &&
          identical(normalized$gwas, parsed$gwas), "programmatic and YAML normalization agree")

  user_config$genome$chains <- list(hg19_to_hg38 = NULL, hg38_to_hg19 = NULL)
  yaml::write_yaml(user_config, yaml_file)
  check(identical(env$read_EPIC2_config(yaml_file, check_files = TRUE)$genome$build, "hg38"),
        "file checking needs no chain for empty legacy mappings")
  user_config$genome$chains <- list(hg19_to_hg38 = "nonexistent.chain.gz")
  yaml::write_yaml(user_config, yaml_file)
  expect_error(env$read_EPIC2_config(yaml_file, check_files = TRUE),
               c("no automatic", "remove genome.chains"))

  # Prove missing required build is rejected before defaults are consulted.
  saved_defaults <- env$EPIC2_defaults
  env$EPIC2_defaults <- function() stop("DEFAULTS_MUST_NOT_RUN")
  legacy_yaml <- list(genome = list(target = "hg38"), atac = list(genome = "hg38"),
                      gwas = list(genome = "hg38"))
  yaml::write_yaml(legacy_yaml, yaml_file)
  expect_error(env$read_EPIC2_config(yaml_file), "genome.build")
  yaml::write_yaml(list(genome = list(build = "")), yaml_file)
  expect_error(env$read_EPIC2_config(yaml_file), "Invalid genome.build")
  env$EPIC2_defaults <- saved_defaults

  cat("PASS:", checks, "common-build configuration checks; no pipeline analysis was run.\n")
  invisible(checks)
}

run_common_build_config_tests(test_package_dir)
