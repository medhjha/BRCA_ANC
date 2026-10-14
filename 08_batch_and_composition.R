# Two checks on the subtype-independent edges.
# First, collection site. Because correlation is computed within each ancestry
# group, we can regress site out within EUR and within AFR separately without
# touching the between-group contrast. Second, composition. We add fibroblast and
# CD8 T-cell marker scores as covariates, beyond the ESTIMATE scores already used,
# to separate coordination from cell-type abundance.

library(data.table)
library(yaml)

cfg <- read_yaml("config.yaml")
inputs <- file.path(cfg$results, "r2_inputs")

Xall <- readRDS(file.path(inputs, "tcga_expr_resid.rds"))
meta <- as.data.frame(readRDS(file.path(inputs, "tcga_meta.rds")))
covs <- readRDS(file.path(inputs, "tcga_covariates.rds"))
edges <- fread(file.path(cfg$results, "subtype_independent_edges.csv"))

meta <- meta[match(colnames(Xall), meta$sample), ]
keep <- meta$ancestry %in% c("EUR", "AFR")
Xall <- Xall[, keep]
meta <- meta[keep, ]
lab <- meta$ancestry
pam <- meta$pam50
cv <- covs[match(colnames(Xall), covs$sample), ]
meta$site <- substr(meta$sample, 6, 7)

fisher_z <- function(r) atanh(pmin(pmax(r, -0.999999), 0.999999))
concordance <- function(a, b) mean(sign(a) == sign(b), na.rm = TRUE)
edge_genes <- unique(c(edges$gene1, edges$gene2))
edge_genes <- edge_genes[edge_genes %in% rownames(Xall)]

# residualize within a group on a design, return residuals and k
residualize_on <- function(samples, design) {
  design <- design[, qr(design)$pivot[seq_len(qr(design)$rank)], drop = FALSE]
  list(resid = t(qr.resid(qr(design), t(Xall[edge_genes, samples]))),
       k = ncol(design))
}

# recompute the differential z for each edge under a within-group adjustment
adjusted_z <- function(res_eur, res_afr, df_eur, df_afr) {
  g2r <- rownames(res_eur)
  sapply(seq_len(nrow(edges)), function(i) {
    g1 <- edges$gene1[i]; g2 <- edges$gene2[i]
    if (!(g1 %in% g2r) || !(g2 %in% g2r)) return(NA)
    (fisher_z(cor(res_afr[g1, ], res_afr[g2, ])) -
       fisher_z(cor(res_eur[g1, ], res_eur[g2, ]))) /
      sqrt(1 / df_eur + 1 / df_afr)
  })
}

report <- function(z, label) {
  ok <- is.finite(z) & is.finite(edges$z)
  n_sig <- sum(p.adjust(2 * pnorm(-abs(z[ok])), "BH") < 0.05)
  cat(sprintf("%s: %d of %d still significant (%.0f%%), direction kept %.0f%%\n",
              label, n_sig, sum(ok), 100 * n_sig / sum(ok),
              100 * concordance(z[ok], edges$z[ok])))
}

# check 1: within-group site regression
site_design <- function(samples) {
  d <- meta[match(samples, meta$sample), ]
  counts <- table(d$site)
  d$site2 <- ifelse(d$site %in% names(counts)[counts >= 3], d$site, "other")
  if (length(unique(d$site2)) < 2) return(matrix(1, nrow(d), 1))
  model.matrix(~ factor(site2), data = d)
}
eur <- meta$sample[lab == "EUR"]
afr <- meta$sample[lab == "AFR"]
re <- residualize_on(eur, site_design(eur))
ra <- residualize_on(afr, site_design(afr))
z_site <- adjusted_z(re$resid, ra$resid, length(eur) - re$k - 3, length(afr) - ra$k - 3)
report(z_site, "within-group site regression")

# check 2: add fibroblast and CD8 marker scores beyond ESTIMATE
fib_markers <- intersect(c("COL1A1", "COL1A2", "COL3A1", "FAP", "PDGFRB",
                           "ACTA2", "POSTN", "LUM", "DCN", "THBS2"), rownames(Xall))
cd8_markers <- intersect(c("CD8A", "CD8B", "GZMA", "GZMK", "NKG7", "PRF1"),
                         rownames(Xall))
fib_score <- colMeans(Xall[fib_markers, ])
cd8_score <- colMeans(Xall[cd8_markers, ])

comp_design <- function(samples) {
  idx <- match(samples, colnames(Xall))
  d <- data.frame(pam = droplevels(factor(pam[idx])),
                  purity = cv$purity[idx], adipo = cv$adipocyte[idx],
                  stromal = cv$StromalScore[idx], immune = cv$ImmuneScore[idx],
                  fib = fib_score[idx], cd8 = cd8_score[idx])
  model.matrix(~ pam + purity + adipo + stromal + immune + fib + cd8, data = d)
}
# exclude the marker genes themselves to avoid circularity
keep_edge <- !(edges$gene1 %in% c(fib_markers, cd8_markers) |
               edges$gene2 %in% c(fib_markers, cd8_markers))
re <- residualize_on(eur, comp_design(eur))
ra <- residualize_on(afr, comp_design(afr))
z_comp <- adjusted_z(re$resid, ra$resid, length(eur) - re$k - 3, length(afr) - ra$k - 3)
z_comp[!keep_edge] <- NA
report(z_comp, "fibroblast and CD8 adjustment")
cat("high survival means coordination, not cell-type abundance\n")
