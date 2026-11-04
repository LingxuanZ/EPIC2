#!/usr/bin/env Rscript
# Run from the directory containing EPIC2_vignette.Rmd.
# Knit the Rmd once, running its stage calls; convert that same knitted text to Markdown.
if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("Install rmarkdown first.")
pandoc_dir <- "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64"
if (!nzchar(Sys.getenv("RSTUDIO_PANDOC")) && file.exists(file.path(pandoc_dir, "pandoc"))) {
  Sys.setenv(RSTUDIO_PANDOC = pandoc_dir)
}
rmarkdown::render("EPIC2_vignette.Rmd", output_format = "html_document",
                  output_file = "EPIC2_vignette.html", quiet = FALSE)
template <- system.file("rmarkdown/templates/github_document/resources/default.md",
                        package = "rmarkdown")
converted <- tempfile(fileext = ".md")
rmarkdown::pandoc_convert("EPIC2_vignette.md", from = "markdown", to = "gfm",
                          output = converted,
                          options = c("--standalone", "--template", template,
                                      "--wrap=none", "--citeproc",
                                      "--bibliography=references_vignette.bib",
                                      "--metadata=link-citations:true"))
stopifnot(file.copy(converted, "EPIC2_vignette.md", overwrite = TRUE))
unlink(converted)
