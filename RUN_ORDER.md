# Run order

Run the scripts in order. Each one reads files written by an earlier script and
writes its own, so they can be run one at a time. Set the paths in config.yaml
first.

## Build the data

    python/parse_ancestry.py          parse the Carrot-Zhang ancestry table
    R/01_download_tcga.R              download RNA, methylation, copy number
    R/02_cohort_audit.R               build the three-layer cohort
    R/03_preprocess.R                 VST, M-values, tumor-half copy number
    R/03b_covariates_and_batch.R      PAM50, purity, ESTIMATE, xCell, batch check
    R/04_wgcna_and_preservation.R     network and module preservation
    R/05_differential_coexpression.R  DGCA, the 4,588 differential edges
    R/05b_build_analysis_inputs.R     write the input objects for scripts 06 to 12
    python/cptac_setup.py             pull the CPTAC RNA, protein, and ancestry

## The analyses that produce the final numbers

    R/06_entanglement_null.R          ancestry vs subtype, the 0.25 result
    R/07_df_correction.R              degrees-of-freedom correction, 3,077 to 1,597
    R/08_batch_and_composition.R      within-group site regression; fibroblast/CD8
    R/09_cptac_reproduction.R         reproduction in CPTAC, by direction
    R/10_enrichment.R                 GO on the two directions
    R/11_panel_invariance.R           panel coherence and the CFA fit check
    R/12_protein_corroboration.R      protein check on the top stromal edges

## Figures

    figures/make_main_figures.R       Figures 1 to 4

## Notes

Script 05b writes the objects (expression matrix, metadata, covariates, pair list)
that scripts 06 to 12 read, so run it after 05. The CPTAC scripts need the CPTAC
inputs built by python/cptac_setup.py; see python/README_cptac.md for the one-time
conversion from CSV to rds.

Script 07 must run before 08 to 12, because those read the subtype-independent edge
set it writes.

The permutation counts are set in config.yaml. The subtype-matched null uses 200
permutations (params: perm_subtype_null), which matches the value reported in the
paper.
