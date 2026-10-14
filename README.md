# Ancestry-differential co-expression in breast cancer

Code for the analysis of ancestry-associated differential gene co-expression in
breast cancer, comparing European-ancestry (EUR) and African-ancestry (AFR)
tumors in TCGA, with reproduction tested in CPTAC.

## What the analysis does

We test whether gene co-expression in breast cancer differs between EUR and AFR
patients, and whether any difference is more than a reflection of molecular
subtype. The short version of the result: the large co-expression modules are the
same in both groups. At the level of individual gene pairs, ancestry produces no
more differential edges than subtype alone, so most of the apparent difference is
explained by subtype. After adjusting within each ancestry group for subtype and
composition, and correcting the test statistic for the degrees of freedom this
uses, a smaller set of about 1,597 edges remains. These split into a stromal
program that is stronger in EUR tumors and reproduces in CPTAC, and an immune
program that is stronger in AFR tumors and does not reproduce.

## Layout

    R/          R scripts for each stage, numbered in run order
    python/     Python scripts (ancestry parsing, CPTAC data)
    figures/    R scripts that make the figures
    config.yaml example config; every script reads paths and parameters from here
    ARCHITECTURE.md  the full analysis flow as a diagram

## How it is organized

The pipeline runs in stages. Each script reads its inputs from disk (saved by an
earlier stage) and writes its outputs back to disk, so stages can be run one at a
time. Paths and parameters live in config.yaml; change the paths there to match
where your data lives. The full run order is in RUN_ORDER.md.

    01 to 03    download TCGA data, build and freeze the cohort, preprocess
    03b         covariates (PAM50, purity, ESTIMATE, xCell) and the batch check
    04          build the co-expression network (WGCNA) and test module preservation
    05          differential co-expression (DGCA), the discovery edges
    05b         write the input objects used by the analyses below
    06          entanglement test: ancestry against a subtype-matched null
    07          degrees-of-freedom correction, the subtype-independent edge set
    08          within-group site regression and the composition check
    09          reproduction in CPTAC, split by direction
    10          Gene Ontology enrichment of the two directions
    11          signature panel coherence and the factor-analysis fit check
    12          protein corroboration of the top stromal edges

## The main result in one line

Most ancestry-associated differential co-expression reflects molecular subtype; a
smaller stromal component is subtype-independent and reproduces in CPTAC, while an
immune component does not. Comparisons of this kind need subtype-matched null
models and a degrees-of-freedom correction.

## Requirements

R 4.4 with Bioconductor 3.20. Main packages: TCGAbiolinks, WGCNA, DGCA, DESeq2,
limma, sesame, genefu, STRINGdb, lionessR, clusterProfiler, estimate,
immunedeconv, lavaan. Python 3.11 with cptac, pandas, numpy, scipy.

Exact versions are listed in the manuscript methods and in Supplementary Table S1.
Three packages are not on CRAN or Bioconductor: DGCA from GitHub
(andymckenzie/DGCA), estimate from R-Forge, and immunedeconv from GitHub
(omnideconv/immunedeconv).

## Repository URL

Replace the placeholder repository URL in the manuscript methods (software and
data availability) with this repository address once it is published.

## Data

All input data is public. TCGA breast cancer data comes from the NCI Genomic Data
Commons through TCGAbiolinks. CPTAC breast cancer data comes from the cptac Python
package. Ancestry calls come from Carrot-Zhang et al. 2020.

## A note on the degrees-of-freedom correction

Residualizing expression on covariates uses up degrees of freedom. When the
differential-correlation test treats the residuals as if the full sample size
remained, the test statistic is too confident, more so in the smaller AFR group.
Scripts 07 onward set the Fisher z variance to 1 over (n minus k minus 3), where k
is the number of covariates. This matters: without it the subtype-independent
count is 3,077, with it the count is 1,597.

## A note on the filesystem

The scripts were written to run off an external drive, so paths point at a mount
point set in config.yaml. On some filesystems (for example exFAT) conda
environments and R package libraries do not work reliably, so keep those on a
local disk and keep only data and results on the external drive if you use one.
