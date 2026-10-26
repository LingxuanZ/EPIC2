% Maintained against the exported function signature in R/genotype.R.
\name{normalize_genotype}
\alias{normalize_genotype}
\title{Normalize common genotype encodings to an allele-dosage matrix}
\description{
Normalize common genotype encodings to an allele-dosage matrix.
}
\usage{
normalize_genotype(genotype, variant_ids = NULL, format = "auto", ploidy = 2, 
    orientation = "auto", missing_values = -9, variant_id_column = NULL) 
}
\arguments{
  \item{genotype}{Matrix-like genotype input, a three-dimensional probability array, or a list of probability-state matrices.}
  \item{variant_ids}{Optional GWAS variant IDs used only to infer orientation. They never supply missing genotype identifiers. Explicit reference-variant identifiers must be present in input row or column names, array dimnames, or the specified variant ID column.}
  \item{format}{One of \code{auto}, \code{dosage}, \code{hard_call}, \code{gt}, \code{allele}, \code{gp}, \code{gl}, \code{pl}, or \code{carrier}.}
  \item{ploidy}{Positive integer copy number.}
  \item{orientation}{Whether variants are rows, columns, or inferred from IDs.}
  \item{missing_values}{Numeric or character missing-value codes.}
  \item{variant_id_column}{Optional variant ID column for data-frame input.}
}
\details{
The result represents allele dosage, not carrier status; carrier-only input is rejected because it cannot recover allele frequency. Genotype encodings must use the allele orientation expected by the GWAS. This function does not infer strand or effect-allele harmonization.
}
\value{
A numeric variant-by-sample allele-dosage matrix, with variant and sample dimnames and a \code{genotype_format} attribute.
}
