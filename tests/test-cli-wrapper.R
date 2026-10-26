# Exercise the repository convenience wrapper in a fresh R process. No analysis
# is run: EPIC2 must reject the missing configuration before entering a stage.
source_root <- Sys.getenv("EPIC2_SOURCE_ROOT", "")
candidates <- c(
  if (nzchar(source_root)) file.path(source_root, "scripts", "run_pipeline.R"),
  file.path("scripts", "run_pipeline.R"),
  file.path("..", "scripts", "run_pipeline.R")
)
wrapper <- candidates[file.exists(candidates)][1L]
if (is.na(wrapper)) {
  message("Repository CLI wrapper test skipped: set EPIC2_SOURCE_ROOT to the source tree.")
} else {
  wrapper <- normalizePath(wrapper, mustWork = TRUE)
  # Carry the library search path into the child, including an isolated test
  # installation supplied by R CMD check or by the developer running this file.
  launcher <- tempfile(fileext = ".R")
  writeLines(c(
    paste0(".libPaths(", paste(capture.output(dput(.libPaths())), collapse = "\n"), ")"),
    paste0("source(", encodeString(wrapper, quote = '"'), ")")
  ), launcher)
  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"), shQuote(launcher), stdout = TRUE, stderr = TRUE
  ))
  unlink(launcher)
  stopifnot(
    identical(attr(output, "status"), 1L),
    any(grepl("A configuration path is required", output, fixed = TRUE)),
    !any(grepl("LDIAG", output, fixed = TRUE))
  )
  message("Repository CLI wrapper reaches EPIC2's configuration validation.")
}
