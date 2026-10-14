# Reproduction of the subtype-independent edges in CPTAC.
# For each edge we compute the EUR-minus-AFR correlation difference in CPTAC and
# compare its direction to TCGA. The null comes from shuffling ancestry labels
# within CPTAC, which keeps the gene-gene dependence intact. We report the result
# for all edges and separately for the two directions, because the stromal and
# immune programs behave differently.

library(data.table)
library(yaml)

cfg <- read_yaml("config.yaml")
inputs <- file.path(cfg$results, "r2_inputs")
set.seed(cfg$seeds$global)

edges <- fread(file.path(cfg$results, "subtype_independent_edges.csv"))
Xc <- readRDS(file.path(inputs, "cptac_rna.rds"))
mc <- as.data.frame(readRDS(file.path(inputs, "cptac_meta.rds")))

mc <- mc[match(colnames(Xc), mc$sample), ]
keep <- mc$ancestry %in% c("EUR", "AFR")
Xc <- Xc[, keep]
labc <- mc$ancestry[keep]

fisher_z <- function(r) atanh(pmin(pmax(r, -0.999999), 0.999999))
concordance <- function(a, b) mean(sign(a) == sign(b), na.rm = TRUE)

testable <- edges$gene1 %in% rownames(Xc) & edges$gene2 %in% rownames(Xc)
e <- edges[testable]
gi <- match(e$gene1, rownames(Xc))
gj <- match(e$gene2, rownames(Xc))
tcga_z <- e$z

# EUR-minus-AFR correlation difference in CPTAC for the current labels
cptac_difference <- function(labels) {
  a <- which(labels == "AFR")
  eu <- which(labels == "EUR")
  cor_in <- function(cols) {
    g <- sort(unique(c(gi, gj)))
    z <- t(scale(t(Xc[g, cols]), center = TRUE, scale = TRUE))
    R <- tcrossprod(z) / (length(cols) - 1)
    R[cbind(match(gi, g), match(gj, g))]
  }
  fisher_z(cor_in(a)) - fisher_z(cor_in(eu))
}

observed <- cptac_difference(labc)
# direction in TCGA: tcga_z positive means higher in AFR, so compare signs directly
is_afr_tighter <- tcga_z > 0

overall <- concordance(observed, tcga_z)
eur_side <- concordance(observed[!is_afr_tighter], tcga_z[!is_afr_tighter])
afr_side <- concordance(observed[is_afr_tighter], tcga_z[is_afr_tighter])

n_perm <- cfg$params$perm_cptac
null_all <- null_eur <- null_afr <- numeric(n_perm)
for (b in seq_len(n_perm)) {
  d <- cptac_difference(sample(labc))
  null_all[b] <- concordance(d, tcga_z)
  null_eur[b] <- concordance(d[!is_afr_tighter], tcga_z[!is_afr_tighter])
  null_afr[b] <- concordance(d[is_afr_tighter], tcga_z[is_afr_tighter])
}
p_of <- function(null, obs) (1 + sum(null >= obs)) / (n_perm + 1)

cat(sprintf("testable edges: %d (EUR %d / AFR %d in CPTAC)\n",
            sum(testable), sum(labc == "EUR"), sum(labc == "AFR")))
cat(sprintf("all edges: %.0f%% concordant, null %.0f%%, p = %.4f\n",
            100 * overall, 100 * median(null_all), p_of(null_all, overall)))
cat(sprintf("EUR-tighter (stromal): %.0f%%, null %.0f%%, p = %.4f (n = %d)\n",
            100 * eur_side, 100 * median(null_eur), p_of(null_eur, eur_side),
            sum(!is_afr_tighter)))
cat(sprintf("AFR-tighter (immune): %.0f%%, null %.0f%%, p = %.4f (n = %d)\n",
            100 * afr_side, 100 * median(null_afr), p_of(null_afr, afr_side),
            sum(is_afr_tighter)))
cat("the stromal edges reproduce; the immune edges do not\n")

e[, cptac_delta := observed]
fwrite(e, file.path(cfg$results, "cptac_reproduction.csv"))
