# Main figures.
# Figure 1: module preservation. Figure 2: the volcano of differential edges.
# Figure 3: reproduction of the subtype-independent edges in CPTAC, by direction.

library(ggplot2)
library(data.table)
library(yaml)

cfg <- read_yaml("config.yaml")
figdir <- file.path(cfg$results, "figures")
dir.create(figdir, showWarnings = FALSE)

plain_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9.5, color = "grey40"),
        plot.caption = element_text(size = 8, color = "grey50", hjust = 0))

save_both <- function(plot, name, w, h) {
  ggsave(file.path(figdir, paste0(name, ".png")), plot, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(figdir, paste0(name, ".pdf")), plot, width = w, height = h)
}

# Figure 1: module preservation
preservation <- fread(file.path(cfg$results, "module_preservation_table.csv"))
f1 <- ggplot(preservation, aes(moduleSize, Zsummary)) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "darkgreen") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "orange") +
  geom_point(aes(color = factor(module)), size = 4) +
  labs(title = "Co-expression modules are preserved from EUR to AFR",
       x = "module size (genes)", y = "preservation Zsummary",
       caption = "Modules above Zsummary 10 are strongly preserved.") +
  plain_theme + theme(legend.position = "none")
save_both(f1, "module_preservation", 7, 5.5)

# Figure 2: volcano, coloured by which group has the tighter correlation
dc <- as.data.table(readRDS(file.path(cfg$results, "dgca_final_EUR_vs_AFR.rds")))
dc[, neglog_fdr := -log10(pmax(fdr, 1e-300))]
dc[, group := fifelse(fdr >= 0.05, "not significant",
              fifelse(zScoreDiff > 0, "higher in AFR", "higher in EUR"))]
xmax <- max(abs(dc$zScoreDiff), na.rm = TRUE)
ymax <- max(dc$neglog_fdr, na.rm = TRUE)
f2 <- ggplot(dc, aes(zScoreDiff, neglog_fdr, color = group)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey45") +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
  scale_color_manual(values = c("higher in EUR" = "#2e6da4",
                                "higher in AFR" = "#c0392b",
                                "not significant" = "grey80"),
                     breaks = c("higher in EUR", "higher in AFR")) +
  annotate("text", x = -xmax * 0.5, y = ymax * 0.97, label = "tighter in EUR",
           color = "#2e6da4", size = 3, hjust = 0.5) +
  annotate("text", x = xmax * 0.5, y = ymax * 0.97, label = "tighter in AFR",
           color = "#c0392b", size = 3, hjust = 0.5) +
  labs(title = "Ancestry-differential co-expression edges",
       subtitle = "colour shows which group has the stronger correlation",
       x = "differential correlation z score (EUR vs AFR)",
       y = "-log10 false discovery rate", color = NULL,
       caption = "Dashed line, FDR 0.05.") +
  plain_theme + theme(legend.position = "top")
save_both(f2, "dgca_volcano", 8, 6)

# Figure 3: CPTAC reproduction, coloured by direction
rep <- fread(file.path(cfg$results, "cptac_reproduction.csv"))
rep[, tcga_delta := -z]    # flip so positive means higher in EUR on both axes
rep[, direction := fifelse(z < 0, "EUR-tighter (stromal)", "AFR-tighter (immune)")]
f3 <- ggplot(rep, aes(tcga_delta, cptac_delta, color = direction)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "grey90", alpha = 0.4) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "grey90", alpha = 0.4) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(size = 0.7, alpha = 0.5) +
  scale_color_manual(values = c("EUR-tighter (stromal)" = "#2e6da4",
                                "AFR-tighter (immune)" = "#c0392b")) +
  labs(title = "Reproduction of edge directions in CPTAC",
       subtitle = "stromal edges reproduce; immune edges do not",
       x = "TCGA effect size (EUR minus AFR)",
       y = "CPTAC effect size (EUR minus AFR)", color = NULL,
       caption = "Grey quadrants are concordant.") +
  plain_theme + theme(legend.position = "top")
save_both(f3, "cptac_reproduction", 7, 6)

cat("main figures written to", figdir, "\n")
