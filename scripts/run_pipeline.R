#!/usr/bin/env Rscript

if (!requireNamespace("EPIC2", quietly = TRUE)) {
  stop("Install the package first with: R CMD INSTALL .", call. = FALSE)
}
LDIAG::ldiag_cli(c("run", commandArgs(trailingOnly = TRUE)))
