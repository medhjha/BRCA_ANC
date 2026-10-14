# Protein corroboration of the top stromal edges.
# For the top 50 EUR-tighter edges with both proteins measured in CPTAC, we check
# whether the direction of the EUR-minus-AFR correlation difference at the protein
# level matches the direction at the RNA level in TCGA. With only 20 AFR samples
# this is a supporting check, not a confirmation.

library(data.table)
library(yaml)

cfg <- read_yaml("config.yaml")
inputs <- file.path(cfg$results, "r2_inputs")
set.seed(cfg$seeds$global)

edges <- fread(file.path(cfg$results, "subtype_independent_edges.csv"))
protein_file <- file.path(inputs, "cptac_protein.rds")
if (!file.exists(protein_file)) {
  stop("cptac_protein.rds not found. Run python/cptac_setup.py with proteomics first.")
}
Xp <- readRDS(protein_file)
mp <- as.data.frame(readRDS(file.path(inputs, "cptac_meta.rds")))
mp <- mp[match(colnames(Xp), mp$sample), ]
keep <- mp$ancestry %in% c("EUR", "AFR")
Xp <- Xp[, keep]
labp <- mp$ancestry[keep]

fisher_z <- function(r) atanh(pmin(pmax(r, -0.999999), 0.999999))
concordance <- function(a, b) mean(sign(a) == sign(b), na.rm = TRUE)

# the top 50 EUR-tighter edges (z < 0 is tighter in EUR) with both proteins present
eur_edges <- edges[z < 0][order(z)]
eur_edges <- eur_edges[gene1 %in% rownames(Xp) & gene2 %in% rownames(Xp)]
top <- eur_edges[1:min(50, nrow(eur_edges))]
cat("top EUR-tighter edges with both proteins:", nrow(top), "\n")

protein_delta <- function(labels) {
  a <- which(labels == "AFR")
  e <- which(labels == "EUR")
  sapply(seq_len(nrow(top)), function(i) {
    fisher_z(cor(Xp[top$gene1[i], a], Xp[top$gene2[i], a])) -
      fisher_z(cor(Xp[top$gene1[i], e], Xp[top$gene2[i], e]))
  })
}

tcga_direction <- top$z
observed <- concordance(protein_delta(labp), tcga_direction)

n_perm <- cfg$params$perm_cptac
null <- replicate(n_perm, concordance(protein_delta(sample(labp)), tcga_direction))
p_value <- (1 + sum(null >= observed)) / (n_perm + 1)

cat(sprintf("protein concordance with TCGA direction: %.0f%%\n", 100 * observed))
cat(sprintf("permutation null: %.0f%% [%.0f to %.0f], p = %.3f\n",
            100 * median(null), 100 * quantile(null, 0.025),
            100 * quantile(null, 0.975), p_value))
cat("this is a supporting check; with 20 AFR samples it is underpowered\n")
