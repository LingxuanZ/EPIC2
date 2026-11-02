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
  \item{missing_values}{Numeric or character missing-value codes. For \code{gl}, numeric codes apply only to entire scalar text genotype cells, never to individual likelihood states; use \code{NA} for missing array/list entries.}
  \item{variant_id_column}{Optional variant ID column for data-frame input.}
}
\details{
The result represents allele dosage, not carrier status; carrier-only input is rejected because it cannot recover allele frequency. Genotype encodings must use the allele orientation expected by the GWAS. This function does not infer strand or effect-allele harmonization.

GL entries are log10 likelihoods, so finite negative values such as \code{-9} are valid states. For example, \code{"0,-9,-18"} is retained with the default missing-value settings. A whole text cell \code{"-9"} is still treated as missing; a three-state vector \code{c(-9, -9, -9)} is an uninformative, valid likelihood vector. Numeric missing codes are not applied to individual GL states, even when explicitly configured. Represent missing numeric array/list genotypes with \code{NA}; character arrays also accept the standard nonnumeric missing tokens and configured nonnumeric codes. For dosage, GP, and PL, configured missing codes apply to numeric entries.
}
\value{
A numeric variant-by-sample allele-dosage matrix, with variant and sample dimnames and a \code{genotype_format} attribute. Allele-label input additionally records \code{counted_alleles} and \code{observed_alleles}; \code{harmonize_gwas()} uses these to orient dosage to the GWAS effect allele. Numerical, GT and probability inputs must already use the GWAS effect-allele orientation.
}
