# Analysis architecture

This is the full flow of the analysis, from raw data to the findings in the paper.
Each stage reads files written by an earlier stage and writes its own outputs, so
stages can be run one at a time.

```mermaid
flowchart TD
    A[TCGA-BRCA download<br/>RNA, methylation, copy number<br/>+ Carrot-Zhang ancestry] --> B[Cohort audit<br/>539 EUR / 110 AFR / 45 admixed]
    B --> C[Preprocess<br/>VST expression, M-values,<br/>tumor-half copy number]
    C --> D[Covariates<br/>PAM50, purity, ESTIMATE, xCell]
    D --> E[Batch check<br/>site/plate vs ancestry<br/>Cramer V 0.58, no ComBat]

    C --> F[WGCNA on EUR<br/>7 modules]
    F --> G[Module preservation EUR to AFR<br/>all preserved, Zsummary 23 to 60]

    C --> H[Clean gene universe<br/>remove low-expression artifact]
    H --> I[DGCA on 539,054 pairs<br/>permutation FDR<br/>4,588 differential edges]

    I --> J[Entanglement test<br/>shuffle ancestry within subtype<br/>5,321 vs null 4,443, p=0.25]

    I --> K[Residualize within group<br/>subtype + composition<br/>degrees-of-freedom correction]
    K --> L[Subtype-independent edges<br/>1,597 df-corrected<br/>vs null 395, p=0.02]

    L --> M[Within-group site regression<br/>75% survive, 99% direction]
    L --> N[Composition check<br/>add fibroblast + CD8<br/>97% survive]
    L --> O[Directional GO enrichment<br/>EUR-tighter: matrix/collagen<br/>AFR-tighter: T-cell/immune]

    L --> P[CPTAC reproduction<br/>permutation null]
    P --> Q[EUR-tighter stromal: 68%, p=0.015<br/>AFR-tighter immune: 51%, p=0.61]

    O --> R[Signature panel invariance<br/>T-cell-inflamed GEP + matrisome<br/>no coherence difference detected]

    classDef neg fill:#fdf1e0,stroke:#c08a2e
    classDef pos fill:#dfe7ef,stroke:#2e4a6b
    classDef key fill:#fbeaea,stroke:#b03030
    class G neg
    class J,Q key
    class L,P pos
```

## Reading the diagram

The two red boxes are the results that shape the paper. The entanglement test
shows that ancestry does not produce more differential edges than subtype alone.
The CPTAC reproduction shows that the stromal edges hold up in an independent
cohort while the immune edges do not.

The blue boxes are the subtype-independent set and its reproduction. The amber
box is the module-level negative result that frames the rest of the work.
