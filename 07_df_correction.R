# The degrees-of-freedom correction and the subtype-independent edge set.
# We residualize expression within each ancestry group on subtype and composition,
# then test all pairs for differential correlation. The test statistic is corrected
# for the number of covariates used, because residualizing uses degrees of freedom
# and an uncorrected test is too confident. We then hold the resulting count to a
# subtype-matched null, the same standard as the raw count.

library(data.table)
library(yaml)

cfg <- read_yaml("config.yaml")
inputs <- file.path(cfg$results, "r2_inputs")
set.seed(cfg$seeds$global)

Xall <- readRDS(file.path(inputs, "tcga_expr_resid.rds"))
meta <- as.data.frame(readRDS(file.path(inputs, "tcga_meta.rds")))
covs <- readRDS(file.path(inputs, "tcga_covariates.rds"))
pairs <- readRDS(file.path(inputs, "tested_pairs.rds"))

meta <- meta[match(colnames(Xall), meta$sample), ]
keep <- meta$ancestry %in% c("EUR", "AFR")
Xall <- Xall[, keep]
meta <- meta[keep, ]
lab <- meta$ancestry
pam <- meta$pam50
cv <- covs[match(colnames(Xall), covs$sample), ]

fisher_z <- function(r) atanh(pmin(pmax(r, -0.999999), 0.999999))

covar_frame <- data.frame(
  pam50 = pam, purity = cv$purity, adipo = cv$adipocyte,
  stromal = cv$StromalScore, immune = cv$ImmuneScore)
complete <- complete.cases(covar_frame) & !is.na(pam)

genes <- unique(c(pairs$gene1, pairs$gene2))
genes <- genes[genes %in% rownames(Xall)]
rows_valid <- which(pairs$gene1 %in% genes & pairs$gene2 %in% genes)

# residualize the edge genes within a set of samples, return residuals and the
# number of covariate columns (k) so the variance can be corrected
residualize <- function(samples) {
  d <- covar_frame[match(samples, colnames(Xall)), ]
  d$pam50 <- droplevels(factor(d$pam50))
  design <- model.matrix(~ pam50 + stromal + immune + purity + adipo, data = d)
  list(resid = t(qr.resid(qr(design), t(Xall[genes, samples]))), k = ncol(design))
}

# differential-correlation z for every pair, corrected for degrees of freedom
count_significant <- function(res_eur, res_afr, df_eur, df_afr) {
  z <- rep(NA_real_, nrow(pairs))
  for (chunk in split(rows_valid, ceiling(seq_along(rows_valid) / 50000))) {
    i <- match(pairs$gene1[chunk], rownames(res_eur))
    j <- match(pairs$gene2[chunk], rownames(res_eur))
    ok <- !is.na(i) & !is.na(j)
    re <- mapply(function(a, b) cor(res_eur[a, ], res_eur[b, ]), i[ok], j[ok])
    ra <- mapply(function(a, b) cor(res_afr[a, ], res_afr[b, ]), i[ok], j[ok])
    z[chunk[ok]] <- (fisher_z(ra) - fisher_z(re)) / sqrt(1 / df_eur + 1 / df_afr)
  }
  list(z = z, n_sig = sum(p.adjust(2 * pnorm(-abs(z)), "BH") < 0.05, na.rm = TRUE))
}

eur_samples <- intersect(meta$sample[lab == "EUR"], colnames(Xall)[complete])
afr_samples <- intersect(meta$sample[lab == "AFR"], colnames(Xall)[complete])
RE <- residualize(eur_samples)
RA <- residualize(afr_samples)
n_eur <- length(eur_samples)
n_afr <- length(afr_samples)

# naive count, then corrected count
naive <- count_significant(RE$resid, RA$resid, n_eur - 3, n_afr - 3)
corrected <- count_significant(RE$resid, RA$resid, n_eur - RE$k - 3, n_afr - RA$k - 3)
cat("uncorrected count (n - 3):", naive$n_sig, "\n")
cat("df-corrected count (n - k - 3):", corrected$n_sig, "\n")

# save the corrected subtype-independent edge set
fdr <- p.adjust(2 * pnorm(-abs(corrected$z)), "BH")
sig <- which(fdr < 0.05 & is.finite(corrected$z))
edges <- data.table(gene1 = pairs$gene1[sig], gene2 = pairs$gene2[sig],
                    z = corrected$z[sig], fdr = fdr[sig])
fwrite(edges, file.path(cfg$results, "subtype_independent_edges.csv"))
cat("saved", nrow(edges), "subtype-independent edges\n")

# hold the corrected count to a subtype-matched null
residualize_subset <- function(samples) {
  d <- covar_frame[match(samples, colnames(Xall)), ]
  d$pam50 <- droplevels(factor(d$pam50))
  design <- model.matrix(~ pam50 + stromal + immune + purity + adipo, data = d)
  list(resid = t(qr.resid(qr(design), t(Xall[genes, samples]))), k = ncol(design))
}
n_perm <- cfg$params$perm_subtype_null
observed <- corrected$n_sig
null_counts <- numeric(n_perm)
idx_complete <- which(complete)
for (b in seq_len(n_perm)) {
  shuffled <- lab
  for (st in unique(pam[idx_complete])) {
    g <- idx_complete[pam[idx_complete] == st]
    shuffled[g] <- sample(lab[g])
  }
  pe <- colnames(Xall)[shuffled == "EUR" & complete]
  pa <- colnames(Xall)[shuffled == "AFR" & complete]
  re <- residualize_subset(pe)
  ra <- residualize_subset(pa)
  null_counts[b] <- count_significant(re$resid, ra$resid,
                                      length(pe) - re$k - 3,
                                      length(pa) - ra$k - 3)$n_sig
  if (b %% 50 == 0) cat("  null permutation", b, "of", n_perm,
                        "running median", round(median(null_counts[1:b])), "\n")
}
p_value <- (1 + sum(null_counts >= observed)) / (n_perm + 1)
cat("\nsubtype-matched null median:", median(null_counts), "\n")
cat("observed:", observed, "| p value:", round(p_value, 4), "\n")

saveRDS(list(observed = observed, null = null_counts, p = p_value,
             naive_count = naive$n_sig),
        file.path(cfg$results, "subtype_independent_null.rds"))
