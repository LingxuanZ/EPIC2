available_stages <- function() c("atac", "gwas", "accessibility", "ld", "wrs", "model", "plot")

write_run_manifest <- function(config, stages, traits) {
  ensure_dir(config$output_dir)
  manifest <- c(
    paste0("software: EPIC2 ", utils::packageVersion("EPIC2")),
    paste0("started: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("R: ", R.version.string),
    paste0("platform: ", R.version$platform),
    paste0("stages: ", paste(stages, collapse = ",")),
    paste0("traits: ", paste(traits, collapse = ",")),
    paste0("genome_build: ", config$genome$build),
    "genome_operation: validation_only",
    "conversion_performed: false",
    paste0("seed: ", config$parameters$seed),
    paste0("workers: ", config$parameters$workers),
    paste0("config: ", config$.config_file %||% "<in-memory configuration>"),
    "",
    "sessionInfo:",
    utils::capture.output(utils::sessionInfo())
  )
  writeLines(manifest, file.path(config$output_dir, "run_manifest.txt"))
  effective <- config
  effective$project$output_dir <- ensure_dir(config$output_dir)
  effective$atac$input <- resolve_config_path(config, config$atac$input)
  effective$atac$fragment_files <- lapply(config$atac$fragment_files, function(x) resolve_config_path(config, x))
  for (trait in names(config$gwas$traits)) {
    for (key in c("summary", "genotype")) {
      effective$gwas$traits[[trait]][[key]] <- resolve_config_path(config, config$gwas$traits[[trait]][[key]])
    }
  }
  effective$plots$sldsc_summary <- resolve_config_path(config, config$plots$sldsc_summary)
  effective[c(".config_file", ".config_dir", "output_dir")] <- NULL
  yaml::write_yaml(effective, file.path(config$output_dir, "config.used.yml"))
  save_rds_atomic(config, file.path(config$output_dir, "config.used.rds"))
  history <- ensure_dir(file.path(config$output_dir, "run_history"))
  id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%OS6", tz = "UTC"), "-", Sys.getpid())
  writeLines(manifest, file.path(history, paste0(id, ".txt")))
  save_rds_atomic(config, file.path(history, paste0(id, ".config.rds")))
  invisible(TRUE)
}

#' Run the EPIC2 analysis and publication-figure workflow
#'
#' @param config YAML path or parsed configuration list.
#' @param stages Ordered stage subset. Defaults to the full workflow.
#' @param traits Optional GWAS-name subset retained for backward compatibility.
#' @param gwas_name Optional name or names under `gwas.traits` to process.
#' @return Named list of stage return values.
#' @export
run_EPIC2 <- function(config, stages = available_stages(), traits = NULL, gwas_name = NULL) {
  if (is.character(config) && length(config) == 1L) config <- read_EPIC2_config(config)
  config <- validate_EPIC2_config(config)
  valid <- available_stages()
  unknown <- setdiff(stages, valid)
  if (length(unknown)) stop("Unknown stage(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  stages <- valid[valid %in% stages]
  traits <- resolve_gwas_names(config, traits = traits, gwas_name = gwas_name)
  ensure_dir(config$output_dir)
  write_run_manifest(config, stages, traits)

  runners <- list(
    atac = function() run_atac_stage(config),
    gwas = function() run_gwas_stage(config, traits = traits),
    accessibility = function() run_accessibility_stage(config, traits),
    ld = function() run_ld_stage(config, traits),
    wrs = function() run_wrs_stage(config, traits),
    model = function() run_model_stage(config, traits),
    plot = function() run_plot_stage(config, traits)
  )
  result <- list()
  for (stage in stages) {
    message("\n[EPIC2] Starting stage: ", stage)
    started <- Sys.time()
    result[[stage]] <- runners[[stage]]()
    message(
      "[EPIC2] Finished stage: ", stage, " in ",
      round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2), " minutes"
    )
  }
  invisible(result)
}

cli_usage <- function() {
  paste(
    "EPIC2: LD-aware scATAC-seq/GWAS integration",
    "",
    "Usage:",
    "  EPIC2 validate CONFIG.yml [--check-files]",
    "  EPIC2 run CONFIG.yml [--stages=atac,gwas,accessibility,ld,wrs,model,plot] [--gwas=HDL,LDL]",
    "  EPIC2 stage CONFIG.yml STAGE [--gwas=HDL,LDL]",
    "",
    "Stages: atac, gwas, accessibility, ld, wrs, model, plot",
    sep = "\n"
  )
}

parse_cli_value <- function(args, option) {
  prefix <- paste0("--", option, "=")
  match <- args[startsWith(args, prefix)]
  if (!length(match)) return(NULL)
  strsplit(sub(prefix, "", match[[1L]], fixed = TRUE), ",", fixed = TRUE)[[1L]]
}

#' Command-line entry point
#'
#' @param args Command-line arguments, excluding the R executable.
#' @return Invisible workflow result; exits with an error for invalid commands.
#' @export
EPIC2_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0L || args[[1L]] %in% c("-h", "--help", "help")) {
    cat(cli_usage(), "\n")
    return(invisible(NULL))
  }
  command <- args[[1L]]
  if (!command %in% c("validate", "run", "stage")) stop(cli_usage(), call. = FALSE)
  if (length(args) < 2L) stop("A configuration path is required.\n\n", cli_usage(), call. = FALSE)
  config_path <- args[[2L]]
  if (command == "validate") {
    config <- read_EPIC2_config(config_path, check_files = "--check-files" %in% args)
    cat("Configuration is valid.\nOutput directory: ", config$output_dir, "\n", sep = "")
    return(invisible(config))
  }

  traits <- parse_cli_value(args, "traits")
  gwas_name <- parse_cli_value(args, "gwas")
  if (!is.null(traits) && !is.null(gwas_name)) {
    stop("Use either --traits or --gwas, not both.", call. = FALSE)
  }
  if (command == "stage") {
    if (length(args) < 3L) stop("stage requires a stage name.\n\n", cli_usage(), call. = FALSE)
    stages <- args[[3L]]
  } else {
    stages <- parse_cli_value(args, "stages") %||% available_stages()
  }
  run_EPIC2(config_path, stages = stages, traits = traits, gwas_name = gwas_name)
}
