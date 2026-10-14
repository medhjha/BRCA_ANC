# The entanglement test.
# The question is whether ancestry produces more differential edges than molecular
# subtype produces on its own. We shuffle ancestry labels within each subtype,
# re-select edges each time, and compare the real count to the shuffled counts.
# Shuffling within subtype keeps the subtype-ancestry link intact, so the null
# answers the right question.

library(data.table)
library(yaml)

cfg <- read_yaml("config.yaml")
inputs <- file.path(cfg$results, "r2_inputs")
set.seed(cfg$seeds$global)

Xall <- readRDS(file.path(inputs, "tcga_expr_resid.rds"))
meta <- as.data.frame(readRDS(file.path(inputs, "tcga_meta.rds")))
pairs <- readRDS(file.path(inputs, "tested_pairs.rds"))

meta <- meta[match(colnames(Xall), meta$sample), ]
keep <- meta$ancestry %in% c("EUR", "AFR")
Xall <- Xall[, keep]
meta <- meta[keep, ]
lab <- meta$ancestry
pam <- meta$pam50

# Fisher z transform and a fast correlation for a set of gene pairs
fisher_z <- function(r) atanh(pmin(pmax(r, -0.999999), 0.999999))

pair_correlation <- function(mat, cols, gi, gj) {
  g <- sort(unique(c(gi, gj)))
  centred <- mat[g, cols, drop = FALSE] - rowMeans(mat[g, cols, drop = FALSE])
  norm <- sqrt(rowSums(centred^2))
  norm[norm == 0] <- NA
  z <- centred / norm
  R <- tcrossprod(z)
  R[cbind(match(gi, g), match(gj, g))]
}

# standardised difference in correlation, positive means higher in AFR
z_difference <- function(mat, labels, gi, gj) {
  a <- which(labels == "AFR")
  e <- which(labels == "EUR")
  ra <- pair_correlation(mat, a, gi, gj)
  re <- pair_correlation(mat, e, gi, gj)
  (fisher_z(ra) - fisher_z(re)) / sqrt(1 / (length(a) - 3) + 1 / (length(e) - 3))
}

# shuffle labels within the levels of a grouping variable
shuffle_within <- function(labels, by) {
  out <- labels
  for (level in unique(by)) {
    idx <- which(by == level)
    if (length(idx) > 1) out[idx] <- labels[idx][sample(length(idx))]
  }
  out
}

genes <- unique(c(pairs$gene1, pairs$gene2))
genes <- genes[genes %in% rownames(Xall)]
gi <- match(pairs$gene1, rownames(Xall))
gj <- match(pairs$gene2, rownames(Xall))
valid <- !is.na(gi) & !is.na(gj)

# the real differential edges and the count at the chosen threshold
z_obs <- z_difference(Xall, lab, gi[valid], gj[valid])
z_cut <- 3.5
observed_count <- sum(abs(z_obs) >= z_cut, na.rm = TRUE)
cat("observed edges at |z| >=", z_cut, ":", observed_count, "\n")

# null: shuffle ancestry within subtype, recount
n_perm <- cfg$params$perm_entanglement
null_counts <- numeric(n_perm)
for (b in seq_len(n_perm)) {
  shuffled <- shuffle_within(lab, pam)
  zb <- z_difference(Xall, shuffled, gi[valid], gj[valid])
  null_counts[b] <- sum(abs(zb) >= z_cut, na.rm = TRUE)
  if (b %% 50 == 0) cat("  permutation", b, "of", n_perm, "\n")
}

p_value <- (1 + sum(null_counts >= observed_count)) / (n_perm + 1)
cat("\nwithin-subtype null median:", median(null_counts), "\n")
cat("empirical p value:", round(p_value, 4), "\n")
cat("ancestry does not clearly exceed subtype if p is above 0.05\n")

saveRDS(list(observed = observed_count, null = null_counts, p = p_value),
        file.path(cfg$results, "entanglement_null.rds"))
