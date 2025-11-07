
# Differential Co-Localization Analysis

<!-- badges: start -->
<!-- badges: end -->

Differential Co-Localization analysis (DiCoLo) is a computational
framework designed for capturing condition-specific co-localized gene
patterns within single-cell RNA sequencing data. The major workflow of
DiCoLo comprises the following three main steps:

- Step1. Compute gene-gene OT distance for each condition
- Step2. Construct differential graph operator
- Step3. Identify differential co-localized genes by spectral analysis
- Downstream tasks
  - Assemble significant DiCoLo genes into co-localized gene modules.

<img src="./figures/DiCoLo_workflow.png" width="100%" height="auto" />

## Example tutorial

Please check the DiCoLo tutorial, available in both
[Rmd](https://github.com/KlugerLab/DiCoLo/blob/main/vignettes/DiCoLo_demo.Rmd)
and
[HTML](https://github.com/KlugerLab/DiCoLo/blob/main/vignettes/DiCoLo_demo.html).
The supporting functions used in the tutorial are available in
[`DiCoLo_functions.R`](https://github.com/KlugerLab/DiCoLo/blob/main/vignettes/DiCoLo_functions.R).
