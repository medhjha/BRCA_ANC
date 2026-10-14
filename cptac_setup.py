# Build the CPTAC inputs used by the reproduction and panel scripts.
# Pulls the CPTAC breast RNA matrix and the consensus ancestry calls, restricts to
# EUR and AFR, and writes them as CSV files the R scripts read.

import os
import cptac
import pandas as pd

RESULTS = os.environ.get("RESULTS", "results")
outdir = os.path.join(RESULTS, "r2_inputs")
os.makedirs(outdir, exist_ok=True)

brca = cptac.Brca()
rna = brca.get_dataframe("transcriptomics", "bcm")
if isinstance(rna.columns, pd.MultiIndex):
    rna.columns = rna.columns.get_level_values(0)
rna = rna.loc[:, ~rna.columns.duplicated()]

# consensus ancestry calls shipped with the cptac package, restricted to breast
ancestry_file = os.path.join(
    os.path.dirname(cptac.__file__),
    "data/harmonized-all_cancers/washu_mssm_consensus_ancestries.tsv.gz")
ancestry = pd.read_csv(ancestry_file, sep="\t", compression="gzip")
ancestry = ancestry[ancestry["cancer_type"] == "BRCA"]
ancestry = ancestry.set_index("case_id")["consensus_pred_ancestry"]

meta = pd.DataFrame({"sample": rna.index})
meta["ancestry"] = meta["sample"].map(lambda s: ancestry.get(s))
meta = meta[meta["ancestry"].isin(["EUR", "AFR"])]
print("CPTAC breast: EUR", (meta.ancestry == "EUR").sum(),
      "AFR", (meta.ancestry == "AFR").sum())

# genes by samples, matching the R convention
rna_out = rna.loc[meta["sample"]].T
rna_out.to_csv(os.path.join(outdir, "cptac_rna.csv"))
meta.to_csv(os.path.join(outdir, "cptac_meta.csv"), index=False)
print("wrote cptac_rna.csv and cptac_meta.csv")
print("convert to rds in R with readRDS-friendly format if needed")


# proteomics, used by R/12_protein_corroboration.R
proteomics = brca.get_dataframe("proteomics", "bcm")
if isinstance(proteomics.columns, pd.MultiIndex):
    proteomics.columns = proteomics.columns.get_level_values(0)
proteomics = proteomics.loc[:, ~proteomics.columns.duplicated()]
prot_meta = pd.DataFrame({"sample": proteomics.index})
prot_meta["ancestry"] = prot_meta["sample"].map(lambda s: ancestry.get(s))
prot_meta = prot_meta[prot_meta["ancestry"].isin(["EUR", "AFR"])]
proteomics.loc[prot_meta["sample"]].T.to_csv(os.path.join(outdir, "cptac_protein.csv"))
print("wrote cptac_protein.csv:",
      "EUR", (prot_meta.ancestry == "EUR").sum(),
      "AFR", (prot_meta.ancestry == "AFR").sum())
