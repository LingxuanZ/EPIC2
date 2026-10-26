\name{EPIC2-package}
\alias{EPIC2-package}
\alias{EPIC2}
\docType{package}
\title{EPIC2: LD-aware scATAC-seq and GWAS integration}
\description{
An R workflow for ATAC quality control, GWAS filtering, SNP- and peak-level
accessibility, local LD modeling, Wilcoxon scoring, LD-whitened association
testing, and aggregated Cauchy combination. An optional plot stage generates
figures from saved analysis outputs.
}
\details{
Input GWAS summary statistics, reference genotypes, ATAC peaks, and fragments
must already use the same declared genome assembly and compatible alleles.
The required \code{genome.build} configuration value has no default; supported
assemblies are hg19/GRCh37 and hg38/GRCh38. The main pipeline validates the
common-build contract and never performs automatic liftOver. The standalone
\code{\link{lift_genomic_ranges}} utility is available for explicit preprocessing.

Use \code{\link{read_EPIC2_config}} and \code{\link{run_EPIC2}} for the staged
workflow, or \code{\link{EPIC2_cli}} for its command-line interface. The package
README and example YAML describe the file layout and complete configuration.
Block permutation and S-LDSC estimation are outside this release.
}
\seealso{
\code{\link{run_EPIC2}}, \code{\link{read_EPIC2_config}},
\url{https://github.com/LingxuanZ/EPIC2}
}
\keyword{package}
