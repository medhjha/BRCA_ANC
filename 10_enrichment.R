# Gene Ontology enrichment of the subtype-independent edges, split by direction.
# Genes in edges tighter in EUR tumors are tested separately from genes in edges
# tighter in AFR tumors, with all genes in tested pairs as the background.

library(data.table)
library(yaml)
library(clusterProfiler)
library(org.Hs.eg.db)

cfg <- read_yaml("config.yaml")

edges <- fread(file.path(cfg$results, "subtype_independent_edges.csv"))
pairs <- readRDS(file.path(cfg$results, "r2_inputs", "tested_pairs.rds"))

# z positive means higher in AFR, so z < 0 is tighter in EUR
eur_edges <- edges[z < 0]
afr_edges <- edges[z > 0]
eur_genes <- unique(c(eur_edges$gene1, eur_edges$gene2))
afr_genes <- unique(c(afr_edges$gene1, afr_edges$gene2))
background <- unique(c(pairs$gene1, pairs$gene2))

to_entrez <- function(symbols) {
  na.omit(mapIds(org.Hs.eg.db, symbols, "ENTREZID", "SYMBOL", multiVals = "first"))
}

run_go <- function(genes, label) {
  result <- tryCatch(
    enrichGO(to_entrez(genes), org.Hs.eg.db,
             universe = as.character(to_entrez(background)),
             ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.10,
             readable = TRUE),
    error = function(e) NULL)
  cat("\n", label, "(", length(genes), "genes):\n")
  if (!is.null(result) && nrow(result) > 0) {
    print(head(as.data.table(result)[, .(Description, qvalue)], 8))
    fwrite(as.data.table(result)[, .(Description, GeneRatio, qvalue)],
           file.path(cfg$results, paste0("enrichment_", label, ".csv")))
  } else {
    cat("  no significant enrichment\n")
  }
  result
}

run_go(eur_genes, "EUR_tighter")
run_go(afr_genes, "AFR_tighter")
cat("\nEUR-tighter edges are extracellular matrix; AFR-tighter edges are immune\n")
