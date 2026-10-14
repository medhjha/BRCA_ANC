# Signature panel invariance.
# We ask whether two published panels measure the same thing in both ancestry
# groups. The primary measure is the mean inter-item correlation, with a bootstrap
# confidence interval for the difference. A single-factor confirmatory factor
# analysis is run as a secondary measure, but we check its fit first, because a
# model that does not fit cannot support an invariance claim.

library(data.table)
library(yaml)
library(lavaan)

cfg <- read_yaml("config.yaml")
inputs <- file.path(cfg$results, "r2_inputs")
set.seed(cfg$seeds$global)

Xall <- readRDS(file.path(inputs, "tcga_expr_resid.rds"))
meta <- as.data.frame(readRDS(file.path(inputs, "tcga_meta.rds")))
meta <- meta[match(colnames(Xall), meta$sample), ]
keep <- meta$ancestry %in% c("EUR", "AFR")
Xall <- Xall[, keep]
meta <- meta[keep, ]
lab <- meta$ancestry
pam <- meta$pam50

# the T-cell inflamed gene expression profile (Ayers et al. 2017)
gep <- intersect(c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9",
                   "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3",
                   "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT"), rownames(Xall))

# the matrix panel: the 18 most variable core-matrisome genes in the pooled cohort,
# a rule chosen before the ancestry comparison and blind to ancestry
core_matrisome <- intersect(c("COL1A1", "COL1A2", "COL2A1", "COL3A1", "COL4A1",
  "COL4A2", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3", "COL11A1", "COL12A1",
  "FN1", "FBN1", "LUM", "DCN", "POSTN", "SPARC", "THBS1", "THBS2", "VCAN", "BGN",
  "FBLN1", "LAMA2", "LAMB1", "NID1", "TNC", "ELN", "MMP2", "MMP14", "TIMP1"),
  rownames(Xall))
ecm <- names(sort(apply(Xall[core_matrisome, ], 1, var), decreasing = TRUE))[1:18]

# residualize panel genes on subtype within a group
residualize_subtype <- function(genes, samples) {
  d <- data.frame(pam = droplevels(factor(pam[match(samples, meta$sample)])))
  design <- model.matrix(~ pam, data = d)
  t(qr.resid(qr(design), t(Xall[genes, samples])))
}
mean_inter_item <- function(mat) {
  cm <- cor(t(mat))
  mean(cm[upper.tri(cm)])
}

# primary measure: difference in mean inter-item correlation, bootstrap interval
panel_coherence <- function(genes, label) {
  eur <- meta$sample[lab == "EUR"]
  afr <- meta$sample[lab == "AFR"]
  me <- residualize_subtype(genes, eur)
  ma <- residualize_subtype(genes, afr)
  observed <- mean_inter_item(me) - mean_inter_item(ma)
  boot <- replicate(2000, {
    be <- sample(ncol(me), replace = TRUE)
    ba <- sample(ncol(ma), replace = TRUE)
    mean_inter_item(me[, be]) - mean_inter_item(ma[, ba])
  })
  ci <- quantile(boot, c(0.025, 0.975))
  cat(sprintf("%s: EUR %.2f, AFR %.2f, difference %.3f [%.3f, %.3f]\n",
              label, mean_inter_item(me), mean_inter_item(ma),
              observed, ci[1], ci[2]))
}

cat("mean inter-item correlation, subtype adjusted:\n")
panel_coherence(gep, "T-cell inflamed profile")
panel_coherence(ecm, "matrix panel")

# secondary measure: confirmatory factor analysis, with a fit check first
single_factor_fit <- function(genes, samples, label) {
  d <- as.data.frame(t(Xall[genes, samples]))
  colnames(d) <- paste0("g", seq_along(genes))
  model <- paste0("F =~ ", paste0("g", seq_along(genes), collapse = " + "))
  fit <- tryCatch(cfa(model, d, std.lv = TRUE), error = function(e) NULL)
  if (is.null(fit)) { cat(label, ": did not converge\n"); return() }
  m <- fitMeasures(fit, c("cfi", "rmsea"))
  cat(sprintf("%s: CFI %.3f, RMSEA %.3f, %s\n", label, m[1], m[2],
              if (m[1] > 0.9) "acceptable fit" else "poor fit, do not interpret"))
}
cat("\nconfirmatory factor analysis fit in each group:\n")
single_factor_fit(gep, meta$sample[lab == "EUR"], "GEP in EUR")
single_factor_fit(gep, meta$sample[lab == "AFR"], "GEP in AFR")
single_factor_fit(ecm, meta$sample[lab == "EUR"], "matrix in EUR")
single_factor_fit(ecm, meta$sample[lab == "AFR"], "matrix in AFR")
cat("\nif fit is poor in both groups the panels are not one-dimensional\n")
cat("and the invariance test is not interpretable\n")
