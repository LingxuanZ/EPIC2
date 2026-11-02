cran <- c(
  "Matrix", "yaml", "data.table", "ggplot2", "igraph", "Seurat", "remotes",
  "knitr", "rmarkdown"
)
missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran, repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
bioconductor <- c(
  "BiocParallel", "GenomeInfoDb", "GenomicRanges",
  "IRanges", "Rsamtools", "S4Vectors", "SNPRelate", "rtracklayer",
  "sparseMatrixStats"
)
missing_bioc <- bioconductor[!vapply(bioconductor, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)

if (!requireNamespace("Signac", quietly = TRUE)) {
  install.packages("Signac", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("presto", quietly = TRUE)) {
  remotes::install_github("immunogenomics/presto@b5df6ee6097eb62522f2e557aa93ac21bab2d05f", upgrade = "never")
}
