# Build the input objects used by scripts 06 to 12.
# Scripts 06 onward work on a small set of saved objects (an expression matrix
# restricted to the tested genes, a metadata frame, a covariate frame, the pair
# list, and the CPTAC matrix). This script writes those objects from the outputs
# of the earlier stages, so the reviewer-response analyses can run on their own.

library(data.table)
library(yaml)
library(SummarizedExperiment)

cfg <- read_yaml("config.yaml")
proc <- file.path(cfg$data, "processed")
interim <- file.path(cfg$data, "interim")
results <- cfg$results
outdir <- file.path(results, "r2_inputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# expression matrix (VST), restricted to the cohort, with gene symbols as rownames
vsd <- readRDS(file.path(proc, "rna_vst.rds"))
rna <- readRDS(file.path(interim, "tcga_rna_se.rds"))
ensembl_to_symbol <- setNames(rowData(rna)$gene_name, rownames(rna))

symbols <- ensembl_to_symbol[rownames(vsd)]
keep <- !is.na(symbols) & symbols != ""
expr <- vsd[keep, ]
symbols <- symbols[keep]

# keep the most variable probe per symbol
ord <- order(apply(expr, 1, var), decreasing = TRUE)
expr <- expr[ord, ]
symbols <- symbols[ord]
expr <- expr[!duplicated(symbols), ]
rownames(expr) <- symbols[!duplicated(symbols)]
saveRDS(expr, file.path(outdir, "tcga_expr_resid.rds"))
cat("expression matrix:", nrow(expr), "genes by", ncol(expr), "samples\n")

# metadata: sample, ancestry, pam50, age, stage, site
coh <- as.data.table(arrow::read_parquet(file.path(interim, "brca_cohort_3layer.parquet")))
common <- readRDS(file.path(proc, "cohort_patients_aligned.rds"))
coh <- coh[patient12 %in% common & cohort %in% c("EUR_primary", "AFR_primary")]
coh[, ancestry := fifelse(cohort == "EUR_primary", "EUR", "AFR")]
cov <- as.data.table(readRDS(file.path(proc, "covariates.rds")))

cd <- as.data.table(as.data.frame(colData(rna)))
cd[, patient12 := substr(colnames(rna), 1, 12)]
cd <- cd[!duplicated(patient12)]
age_col <- grep("age_at_diagnosis|age_at_index", names(cd), value = TRUE, ignore.case = TRUE)[1]
stage_col <- grep("ajcc_pathologic_stage", names(cd), value = TRUE)[1]
cd[, age := { v <- as.numeric(get(age_col)); if (median(v, na.rm = TRUE) > 200) v / 365.25 else v }]
cd[, stage := gsub("Stage |A|B|C", "", as.character(get(stage_col)))]

meta <- merge(coh[, .(sample = patient12, ancestry)],
              cov[, .(patient12, pam50 = pam50_genefu)],
              by.x = "sample", by.y = "patient12", all.x = TRUE)
meta <- merge(meta, cd[, .(patient12, age, stage)],
              by.x = "sample", by.y = "patient12", all.x = TRUE)
meta <- as.data.frame(meta[sample %in% colnames(expr)])
saveRDS(meta, file.path(outdir, "tcga_meta.rds"))
cat("metadata:", nrow(meta), "samples\n")

# covariates: sample, purity, adipocyte, StromalScore, ImmuneScore
# the adipocyte score is recomputed with full xCell if it is not stored
adipo_col <- grep("adipo", names(cov), value = TRUE, ignore.case = TRUE)[1]
if (is.na(adipo_col)) {
  library(xCell)
  tpm <- assay(rna, "tpm_unstrand")
  sym <- rowData(rna)$gene_name
  ok <- !is.na(sym) & sym != ""
  tpm <- tpm[ok, ]; rownames(tpm) <- make.unique(sym[ok])
  colnames(tpm) <- substr(colnames(rna), 1, 12)
  tpm <- tpm[, colnames(tpm) %in% meta$sample]
  xc <- xCellAnalysis(tpm)
  adipocyte <- xc["Adipocytes", ][cov$patient12]
} else {
  adipocyte <- cov[[adipo_col]]
}
covars <- data.frame(sample = cov$patient12, purity = cov$CPE,
                     adipocyte = adipocyte,
                     StromalScore = cov$StromalScore,
                     ImmuneScore = cov$ImmuneScore)
covars <- covars[covars$sample %in% colnames(expr), ]
saveRDS(covars, file.path(outdir, "tcga_covariates.rds"))
cat("covariates:", nrow(covars), "samples\n")

# the pre-specified pair list, with gene symbols
pairset <- as.data.table(readRDS(file.path(proc, "dgca_pairset.rds")))
pairset[, g1 := ensembl_to_symbol[Gene1]][, g2 := ensembl_to_symbol[Gene2]]
pairs <- as.data.frame(pairset[!is.na(g1) & !is.na(g2), .(gene1 = g1, gene2 = g2)])
saveRDS(pairs, file.path(outdir, "tested_pairs.rds"))
cat("tested pairs:", nrow(pairs), "\n")

cat("\ninput objects written to", outdir, "\n")
cat("the CPTAC matrix is built separately by python/cptac_setup.py\n")
