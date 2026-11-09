project_dir <- "/path/to/fetal_tissue"
source_dir <- project_dir
output_dir <- file.path(project_dir, "stage1_reproduced")

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomeInfoDb)
  library(GenomicRanges)
  library(Matrix)
  library(dplyr)
  library(EnsDb.Hsapiens.v75)
  library(SNPRelate)
})

traits <- c("HDL", "LDL", "TC", "TG", "SCZ", "SCZ2", "BIP", "SCZBIP")
tissues <- c(
  "heart", "eye", "intestine", "muscle", "adrenal",
  "liver", "kidney", "lung", "pancreas", "cerebrum"
)

read_object <- function(path, name) {
  e <- new.env(parent = baseenv())
  load(path, envir = e)
  if (!exists(name, envir = e, inherits = FALSE)) {
    stop("Object '", name, "' is missing from ", path)
  }
  e[[name]]
}

one_file <- function(files, tissue) {
  selected <- files[grepl(tissue, basename(files), fixed = TRUE)]
  if (length(selected) != 1L) {
    stop("Provide exactly the original selected file for tissue: ", tissue)
  }
  selected
}

outputs <- file.path(
  output_dir,
  c(
    "atac.qc.withFragments.cellQC.peakQC.rds",
    paste0(traits, ".gwas.pruned.noamb.rds"),
    paste0(traits, ".gwas.pruned.noamb.genotype.rda")
  )
)
if (any(file.exists(outputs))) stop("Output files already exist.")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


######==============================
###### ATAC QC
######==============================

# Upstream ATAC preprocessing.
# Published processed ATAC objects and cell annotations: GSE149683.
# Use the same selected fragment file for each tissue, with its .tbi index.

atac_dir <- file.path(source_dir, "atac_fetal")

rds_files <- list.files(
  file.path(atac_dir, "rds"),
  pattern = "\\.RDS$",
  full.names = TRUE
)
fragment_files <- list.files(
  file.path(atac_dir, "fragments"),
  pattern = "\\.txt\\.gz$",
  full.names = TRUE
)
cells_qc <- read.table(
  file.path(
    atac_dir,
    "GSE149683_Metadata_of_high_quality_cells_after_doublet_filtering.txt.gz"
  ),
  sep = "\t"
)
stopifnot("cell" %in% names(cells_qc))

counts_by_tissue <- list()
metadata_by_tissue <- list()

for (tissue in tissues) {
  rds_path <- one_file(rds_files, tissue)
  fragment_path <- one_file(fragment_files, tissue)
  stopifnot(file.exists(paste0(fragment_path, ".tbi")))
  
  atac <- UpdateSeuratObject(readRDS(rds_path))
  
  if (tissue == "cerebrum") {
    atac$sample_name <- gsub("brain", "cerebrum", atac$sample_name)
  }
  
  # Retain the sample corresponding to the selected fragment file.
  selected <- vapply(
    atac$sample_name,
    function(x) grepl(x, basename(fragment_path)),
    logical(1)
  )
  atac <- atac[, selected]
  atac <- atac[, colnames(atac) %in% cells_qc$cell]
  
  chrom_assay <- CreateChromatinAssay(
    counts = GetAssayData(atac, assay = "peaks", layer = "counts"),
    sep = c("-", "-"),
    genome = "hg19",
    min.cells = 0,
    min.features = 100
  )
  atac <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    meta.data = atac[[]]
  )
  
  # Remove uncertain annotations and cell types with fewer than 300 cells.
  atac <- atac[, !grepl("\\?", atac$cell_type)]
  atac <- atac[, !grepl("Unknown", atac$cell_type)]
  atac$cell_type <- droplevels(atac$cell_type)
  
  small_types <- names(which(table(atac$cell_type) < 300))
  atac <- atac[, !atac$cell_type %in% small_types]
  atac$cell_type <- droplevels(atac$cell_type)
  
  atac$tissue_cell_type <- paste0(atac$tissue, "_", atac$cell_type)
  
  counts_by_tissue[[tissue]] <- GetAssayData(
    atac, assay = "ATAC", layer = "counts"
  )
  metadata_by_tissue[[tissue]] <- data.frame(
    cell_type = atac$cell_type,
    tissue = atac$tissue,
    tissue_cell_type = atac$tissue_cell_type,
    row.names = colnames(atac)
  )
}

# Combine tissues on the shared, identically ordered peak set.
stopifnot(all(vapply(
  counts_by_tissue,
  function(x) identical(rownames(x), rownames(counts_by_tissue[[1]])),
  logical(1)
)))

counts <- do.call(cbind, unname(counts_by_tissue))
metadata <- do.call(rbind, unname(metadata_by_tissue))

stopifnot(
  !anyDuplicated(colnames(counts)),
  identical(colnames(counts), rownames(metadata))
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c("-", "-"),
  genome = "hg19",
  min.cells = 0,
  min.features = 500
)
atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC",
  meta.data = metadata
)

rm(counts_by_tissue, metadata_by_tissue, counts, metadata, chrom_assay)

# Initial count- and feature-based cell filtering.
atac <- subset(
  atac,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < median(atac$nCount_ATAC) + 3 * mad(atac$nCount_ATAC) &
    nFeature_ATAC > 1000 &
    nFeature_ATAC < median(atac$nFeature_ATAC) + 3 * mad(atac$nFeature_ATAC)
)

small_groups <- names(which(table(atac$tissue_cell_type) < 500))
atac <- atac[, !atac$tissue_cell_type %in% small_groups]

# LSI and clustering are required for the following cell-selection step.
atac <- RunTFIDF(atac, verbose = FALSE)
atac <- FindTopFeatures(atac, min.cutoff = "q0", verbose = FALSE)
atac <- RunSVD(atac, verbose = FALSE)
atac <- FindNeighbors(
  atac, reduction = "lsi", dims = 2:30, verbose = FALSE
)
atac <- FindClusters(atac, algorithm = 3, verbose = FALSE)

# Within each cluster, retain cells from its most frequent tissue.
cluster_tissue <- table(atac$seurat_clusters, atac$tissue)
major_tissue <- vapply(
  seq_len(nrow(cluster_tissue)),
  function(i) colnames(cluster_tissue)[which.max(cluster_tissue[i, ])],
  character(1)
)
keep_pairs <- paste0(rownames(cluster_tissue), "_", major_tissue)

atac <- atac[
  ,
  paste0(atac$seurat_clusters, "_", atac$tissue) %in% keep_pairs
]

small_tissues <- names(which(table(atac$tissue) < 3000))
atac <- atac[, !atac$tissue %in% small_tissues]

# Additional QC from Fetal_Tissue.Rmd.
DefaultAssay(atac) <- "ATAC"

frag_list <- lapply(unique(as.character(atac$tissue)), function(tissue) {
  CreateFragmentObject(
    path = one_file(fragment_files, tissue),
    cells = colnames(atac)[atac$tissue == tissue],
    validate.fragments = TRUE,
    verbose = FALSE
  )
})
Fragments(atac[["ATAC"]]) <- frag_list

gene_annot <- GetGRangesFromEnsDb(
  ensdb = EnsDb.Hsapiens.v75,
  verbose = FALSE
)
seqlevels(gene_annot) <- paste0("chr", seqlevels(gene_annot))
genome(gene_annot) <- "hg19"
Annotation(atac[["ATAC"]]) <- gene_annot

atac <- NucleosomeSignal(atac, verbose = FALSE)
atac <- TSSEnrichment(atac, verbose = FALSE)

# Preserve the Rmd's FRIP and blacklist calculations.
# FeatureMatrix counts are summed over reduced intervals. A fragment spanning
# separate intervals can still contribute to more than one interval.
# The original denominator calculation requires identity barcode mapping.
stopifnot(all(vapply(frag_list, function(frag) {
  cells <- GetFragmentData(frag, slot = "cells")
  identical(names(cells), unname(cells))
}, logical(1))))

total_frag_list <- lapply(frag_list, function(frag) {
  CountFragments(
    fragments = GetFragmentData(frag, slot = "path"),
    cells = names(GetFragmentData(frag, slot = "cells")),
    verbose = FALSE
  )
})

total_frag_df <- bind_rows(total_frag_list) %>%
  group_by(CB) %>%
  summarise(
    frequency_count = sum(frequency_count),
    .groups = "drop"
  )

total_fragments <- setNames(
  total_frag_df$frequency_count,
  total_frag_df$CB
)
atac$total_unique_fragments <- total_fragments[colnames(atac)]

peak_union <- reduce(
  StringToGRanges(rownames(atac[["ATAC"]]), sep = c("-", "-")),
  ignore.strand = TRUE
)
peak_union_mat <- FeatureMatrix(
  fragments = frag_list,
  features = peak_union,
  cells = colnames(atac),
  verbose = FALSE
)

peak_counts <- Matrix::colSums(peak_union_mat)
atac$unique_fragments_in_peaks <- peak_counts[colnames(atac)]
atac$FRIP <- atac$unique_fragments_in_peaks / atac$total_unique_fragments
rm(peak_union_mat)

data("blacklist_hg19", package = "Signac")

blacklist_mat <- FeatureMatrix(
  fragments = frag_list,
  features = reduce(blacklist_hg19, ignore.strand = TRUE),
  cells = colnames(atac),
  verbose = FALSE
)
blacklist_counts <- Matrix::colSums(blacklist_mat)
atac$blacklist_fragments <- blacklist_counts[colnames(atac)]
atac$BlacklistRatio <- atac$blacklist_fragments / atac$total_unique_fragments
rm(blacklist_mat)

# Cell QC: preserve the sequential filtering order.
# Calculate the first percentile after filtering on TSS enrichment and FRIP.
atac <- subset(atac, subset = TSS.enrichment >= 2.5)
atac <- subset(atac, subset = FRIP >= 0.3)
atac <- subset(
  atac,
  subset = nCount_ATAC >= quantile(nCount_ATAC, 0.01, na.rm = TRUE)
)
atac <- subset(atac, subset = BlacklistRatio < 0.02)
atac <- subset(atac, subset = nucleosome_signal <= 4)

group_n <- table(atac$tissue_cell_type)
keep_groups <- names(group_n[group_n >= 50])
atac <- subset(
  atac,
  cells = colnames(atac)[
    as.character(atac$tissue_cell_type) %in% keep_groups
  ]
)

# Peak QC: remove zero-count peaks, retain autosomes, restrict peak widths,
# and exclude peaks overlapping the hg19 blacklist.
counts <- GetAssayData(atac, assay = "ATAC", layer = "counts")
atac <- atac[Matrix::rowSums(counts) > 0, ]

atac <- subset(
  atac,
  features = rownames(atac)[
    gsub("-.*", "", rownames(atac)) %in% paste0("chr", 1:22)
  ]
)

atac <- subset(
  atac,
  features = rownames(atac)[
    width(granges(atac)) >= 100 & width(granges(atac)) <= 1000
  ]
)

peak_gr <- StringToGRanges(rownames(atac), sep = c("-", "-"))
atac <- subset(
  atac,
  features = GRangesToString(
    peak_gr[!peak_gr %over% blacklist_hg19]
  )
)

saveRDS(
  atac,
  file.path(output_dir, "atac.qc.withFragments.cellQC.peakQC.rds")
)


######==============================
###### GWAS QC
######==============================
gwas_dir <- file.path(source_dir, "gwas")

aligned_files <- c(
  HDL = "HDL_Aligned.rda",
  LDL = "LDL_Aligned.rda",
  TC = "TC_Aligned.rda",
  TG = "TG_Aligned.rda",
  SCZ = "SCZ_aligned.rda",
  SCZ2 = "SCZ2_aligned.rda",
  BIP = "BIP_aligned.rda",
  SCZBIP = "SCZBIP_aligned.rda"
)

prune_snps <- function(gwas, genotype) {
  gds_path <- tempfile(
    pattern = "gwas-pruning-",
    tmpdir = output_dir,
    fileext = ".gds"
  )
  on.exit(unlink(gds_path), add = TRUE)
  
  snpgdsCreateGeno(
    gds_path,
    genmat = genotype,
    sample.id = colnames(genotype),
    snp.id = rownames(genotype),
    snp.chromosome = gwas$chr,
    snp.position = gwas$pos,
    snp.allele = paste0(gwas$A1, "/", gwas$A2),
    snpfirstdim = TRUE
  )
  
  genofile <- snpgdsOpen(gds_path)
  on.exit(snpgdsClose(genofile), add = TRUE, after = FALSE)
  
  # Preserve the original seed, threshold and other package defaults.
  set.seed(1000)
  unlist(
    snpgdsLDpruning(
      genofile,
      ld.threshold = 0.8,
      verbose = FALSE
    )
  )
}

for (trait in traits) {
  gwas <- read_object(
    file.path(gwas_dir, aligned_files[[trait]]),
    "gwas.align"
  )
  genotype <- read_object(
    file.path(gwas_dir, "genotype.rda"),
    "genotype"
  )
  
  # Preserve the original genomic ordering.
  gwas_gr <- GRanges(
    seqnames = gwas$chr,
    ranges = IRanges(start = gwas$pos, end = gwas$pos + 1)
  )
  seqlevelsStyle(gwas_gr) <- "NCBI"
  
  gwas$chr <- as.character(seqnames(gwas_gr))
  gwas$granges.str <- GRangesToString(gwas_gr, sep = c("-", "-"))
  
  gwas_gr <- sort(sortSeqlevels(gwas_gr))
  gwas <- gwas[
    match(
      GRangesToString(gwas_gr, sep = c("-", "-")),
      gwas$granges.str
    ),
    , drop = FALSE
  ]
  rownames(gwas) <- NULL
  
  # Match the reference panel and remove duplicate rsIDs.
  gwas <- gwas[
    as.character(gwas$rsid) %in% rownames(genotype),
    , drop = FALSE
  ]
  genotype <- genotype[
    match(gwas$rsid, rownames(genotype)),
    , drop = FALSE
  ]
  
  unique_snp <- !duplicated(gwas$rsid)
  genotype <- genotype[unique_snp, , drop = FALSE]
  gwas <- gwas[unique_snp, , drop = FALSE]
  
  # Preserve the historical filter.
  # rowMeans(genotype > 0) is carrier frequency, not dosage-based MAF.
  carrier_fraction <- rowMeans(genotype > 0)
  keep <- gwas$MAF >= 0.05 & carrier_fraction >= 0.05
  stopifnot(!anyNA(keep))
  
  genotype <- genotype[keep, , drop = FALSE]
  gwas <- gwas[keep, , drop = FALSE]
  
  # LD pruning.
  snp_set <- prune_snps(gwas, genotype)
  gwas <- gwas[match(snp_set, gwas$rsid), , drop = FALSE]
  genotype <- genotype[
    match(snp_set, rownames(genotype)),
    , drop = FALSE
  ]
  stopifnot(all(gwas$rsid == rownames(genotype)))
  
  # Fetal_Tissue.Rmd: remove strand-ambiguous variants after pruning.
  ambiguous <- (gwas$A1 == "A" & gwas$A2 == "T") |
    (gwas$A1 == "T" & gwas$A2 == "A") |
    (gwas$A1 == "C" & gwas$A2 == "G") |
    (gwas$A1 == "G" & gwas$A2 == "C")
  stopifnot(!anyNA(ambiguous))
  
  gwas <- gwas[!ambiguous, , drop = FALSE]
  genotype <- genotype[!ambiguous, , drop = FALSE]
  
  stopifnot(
    nrow(gwas) == nrow(genotype),
    all(gwas$rsid == rownames(genotype))
  )
  
  saveRDS(
    gwas,
    file.path(output_dir, paste0(trait, ".gwas.pruned.noamb.rds"))
  )
  save(
    genotype,
    file = file.path(
      output_dir,
      paste0(trait, ".gwas.pruned.noamb.genotype.rda")
    )
  )
  
  rm(gwas, genotype)
  invisible(gc())
}
