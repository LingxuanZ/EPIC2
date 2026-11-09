# EPIC2 fetal tissue pipeline: Stage 1 QC and Stages 2–5

Run **Stage 1 first**, then submit **Stages 2–5**. The configured traits are HDL, LDL, TC, TG, SCZ, SCZ2, BIP and SCZBIP. All input coordinates must use hg19/GRCh37.

## 1. Set paths and prepare inputs

Keep `Stage1_QC_ATAC_GWAS.R`, `config.yml`, `run_stages_2_to_5_EPIC2.R` and `submit_stages_2_to_5_EPIC2.sh` together in the HPC script directory. Replace these placeholders before running:

| File | Setting | What to supply |
| --- | --- | --- |
| `Stage1_QC_ATAC_GWAS.R` | `project_dir = "/path/to/fetal_tissue"` | The directory containing `atac_fetal/` and `gwas/`. |
| `config.yml` | `/path/to/fetal_tissue` in all input paths | The same `project_dir`; retain the `stage1_reproduced/` subdirectory. |
| `config.yml` | `/path/to/epic2_results` in `project.output_dir` | The parent directory for Stage 2–5 results. |
| `submit_stages_2_to_5_EPIC2.sh` | `EPIC2_SLURM_ACCOUNT="YOUR_SLURM_ACCOUNT"` | Your Slurm allocation/account. |
| `submit_stages_2_to_5_EPIC2.sh` | `BASE_LIBRARY="/path/to/R/library"` | Your installed R dependency library. |

Use absolute paths accessible from the compute nodes. The supplied shell script uses `GCC/14.2.0` and `R/4.5.1`; adapt the module commands to your cluster if needed. Dependencies must already be installed. Stage 1 loads Seurat, Signac, GenomeInfoDb, GenomicRanges, Matrix, dplyr, EnsDb.Hsapiens.v75 and SNPRelate. The submission script also requires `yaml`, `remotes` and EPIC2's dependencies.

Stage 1 reads these inputs directly below `/path/to/fetal_tissue/`:

- `atac_fetal/rds/`: published tissue-specific `.RDS` objects from GSE149683.
- `atac_fetal/fragments/`: the selected `.txt.gz` fragment files and their `.tbi` indexes.
- `atac_fetal/GSE149683_Metadata_of_high_quality_cells_after_doublet_filtering.txt.gz`: the high-quality cell annotations.
- `gwas/HDL_Aligned.rda`, `LDL_Aligned.rda`, `TC_Aligned.rda` and `TG_Aligned.rda`.
- `gwas/SCZ_aligned.rda`, `SCZ2_aligned.rda`, `BIP_aligned.rda` and `SCZBIP_aligned.rda`.
- `gwas/genotype.rda`: the reference genotype matrix, stored as `genotype`, with rsIDs in row names and samples in columns. Each aligned GWAS file must contain `gwas.align`.

The ATAC script expects exactly one matching RDS and one selected fragment file for each of heart, eye, intestine, muscle, adrenal, liver, kidney, lung, pancreas and cerebrum. Tissue names must appear in the filenames. The GWAS filenames are case-sensitive.

Stage 1 incorporates upstream preprocessing and the additional QC in `Fetal_Tissue.Rmd`. It starts from published processed ATAC objects and already aligned GWAS/reference files. Code that generates the aligned GWAS files and `genotype.rda` from their original releases is not included. The regenerated Stage 1 inputs have not yet been validated against the manuscript inputs; exact reproduction also depends on the original source files, software versions and upstream cell selection.

## 2. Run Stage 1: ATAC and GWAS QC

Run from the script directory on a compute node with sufficient resources. For the supplied module environment:

```bash
module purge
module load GCC/14.2.0 R/4.5.1
export R_LIBS_USER="/path/to/R/library"
export R_LIBS="${R_LIBS_USER}"
Rscript --vanilla Stage1_QC_ATAC_GWAS.R
```

Use the same real R library path as `BASE_LIBRARY` in the submission script.

- **ATAC:** select samples and annotated cells, apply the upstream clustering-based filter, then perform cell and peak QC.
- **GWAS:** match reference genotypes, remove duplicate rsIDs, apply the original frequency filter, run SNPRelate LD pruning and remove strand-ambiguous variants. The original reference filter uses `rowMeans(genotype > 0)`, which is carrier frequency rather than dosage-based MAF.

Stage 1 writes 17 files under `/path/to/fetal_tissue/stage1_reproduced/`:

- `atac.qc.withFragments.cellQC.peakQC.rds`.
- `<TRAIT>.gwas.pruned.noamb.rds` for each of the eight traits.
- `<TRAIT>.gwas.pruned.noamb.genotype.rda` for each trait, containing `genotype`.

These paths match the inputs in `config.yml`. The ATAC object contains `tissue_cell_type`, `tissue` and attached fragment references. Keep the fragment files and indexes accessible after QC. Stage 1 stops if any target output already exists.

Wait for Stage 1 to finish successfully before submitting Stage 2–5. The submission script below does not launch Stage 1.

## 3. Submit Stages 2–5

From the script directory on the login node, run:

```bash
bash submit_stages_2_to_5_EPIC2.sh
```

The shell script reads `config.yml` from its own directory and checks that all configured inputs are readable. It installs the pinned EPIC2 commit into `~/.local/share/EPIC2/<COMMIT>/R_library` if that copy is not already installed. Dependencies are read from `BASE_LIBRARY` and the system R libraries.

Stage 2 uses the Stage 1 GWAS outputs with `ld_prune: false`, `remove_ambiguous: false`, `maf_min: 0` and `reference_missing_rate_max: 1`. It retains EPIC2's validation of coordinates, identifiers and effect statistics. The accessibility stage removes SNPs without fragment coverage. No cell resampling is performed.

The submission creates one Slurm array with eight tasks. Each trait runs these calls **sequentially**, while different traits can run **in parallel**:

```r
EPIC2::run_gwas_stage(config, traits = trait)
EPIC2::run_accessibility_stage(config, traits = trait)
EPIC2::run_ld_stage(config, traits = trait)
EPIC2::run_wrs_stage(config, traits = trait)
EPIC2::run_model_stage(config, traits = trait)
```

The first two calls prepare Stage 2. A failed call stops that trait before later stages run. Each task requests 4 CPUs, 160 GB RAM and 48 hours. `parameters.workers: 1` keeps accessibility counting within each trait serial; the Slurm array supplies parallelism across traits.

R, shell and YAML files stay in the uploaded directory; the script does not copy them to scratch. Keep the R/YAML files unchanged while jobs are queued or running. Their checksums are recorded and checked when each task starts.

## 4. Stage 2–5 results and logs

The example output root in `config.yml` is:

```text
/path/to/epic2_results/stages_2_to_5/
```

Each `<TRAIT>/` contains:

- `01_atac/atac.qc.rds`: a link to the existing QC object.
- `02_gwas/`: prepared GWAS, reference genotypes and GWAS QC summaries.
- `02_accessibility/`: SNP/peak count matrices and SNP–peak overlaps.
- `03_ld/`: newly calculated SNP covariance, SNP whitening matrix and peak LD model.
- `04_wrs/`: complete SNP/peak WRS tables in RDS format.
- `05_models/`: component and Cauchy results in both RDS and TSV formats, for each annotation level.
- `config.used.rds`, `sessionInfo.txt` and, after all stages succeed, `COMPLETED.txt`.

Slurm logs are in `log/`. Existing submissions and trait outputs are protected. To start another complete run, choose a fresh `project.output_dir` in YAML.

- Each task's `.out` records its trait, Slurm job ID, and UTC `START`/`DONE` messages for setup and every stage. Completed stages report elapsed seconds.
- An R error records `FAILED`, the stage, elapsed time and the original error message. Later stages do not run. The task exit message reports the process exit code.
- EPIC2's original progress messages, warnings and errors remain in `.err`; a nonempty `.err` does not by itself mean failure. Package installation is recorded in `log/installation.log`.
- A forced termination may prevent a final log message. Check the Slurm state and exit code if a task has no final message or `COMPLETED.txt`.

LD uses a 5 Mb window, ridge 0.01 and eigenvalue floor `1e-8`. The dense-block guard is 15,000 features. A larger block stops with an error; the pipeline does not silently split or approximate it. Diagnostic peak checkpoints and WRS TSV exports are disabled to reduce storage; the full model inputs and WRS RDS tables are retained.

Stages 2–5 calculate GWAS/accessibility outputs, LD, WRS and models from the Stage 1 inputs. Plotting and result comparison are separate steps. 
