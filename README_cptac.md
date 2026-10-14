# CPTAC inputs

cptac_setup.py writes three CSV files into results/r2_inputs:
cptac_rna.csv, cptac_meta.csv, and cptac_protein.csv.

The R scripts read .rds versions. Convert them once with this in R:

    outdir <- "results/r2_inputs"
    saveRDS(as.matrix(read.csv(file.path(outdir, "cptac_rna.csv"),
            row.names = 1, check.names = FALSE)),
            file.path(outdir, "cptac_rna.rds"))
    saveRDS(read.csv(file.path(outdir, "cptac_meta.csv")),
            file.path(outdir, "cptac_meta.rds"))
    saveRDS(as.matrix(read.csv(file.path(outdir, "cptac_protein.csv"),
            row.names = 1, check.names = FALSE)),
            file.path(outdir, "cptac_protein.rds"))
