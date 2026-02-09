# N-N-A against H-S

Pipeline for predicting **Half-Sib (H-S)** vs **Not Applicable (N/A)** pairs from genotype/π̂ (pi-hat) features.

## Contents

- **`GTNAall2.R`** — Main R script: data read/merge, filtering, univariate GMMs, supervised classifier (4-D + calibration), multivariate 2-component Gaussian mixture (4-D) → P(H-S), cutoff choice, plots, and output table.

## Usage

```bash
Rscript GTNAall2.R <input_file_path>
```

Example:
```bash
Rscript GTNAall2.R /path/to/pairs_with_type.tsv
```

Input: tab-separated file with columns for pair IDs, `type` (e.g. 100 = N/A, 2/7/11 = H-S), `Pihat`, `Pihat2`, and sum columns `sum_pihat_max_5`, `sum_pihat_remaining_max_5`, `sum_pihat_remaining_min_5`, `sum_pihat_opposite_5`. See script header for full column layout.

## Dependencies

- **R** packages: `mixtools`, `data.table`
- Optional: `glmnet`, `mvtnorm`

Install in R:
```r
install.packages(c("mixtools", "data.table", "mvtnorm"))
# Optional: install.packages("glmnet")
```

## Outputs

The script writes (paths are set inside the script; change them as needed):

- GMM histograms (univariate, Case/Control and All)
- Supervised classifier: calibration, P(HS) by group, ROC
- 4-D mixture: P(HS) histograms, cutoff, calibration, ROC, P(HS) vs age difference
- TSV: P(HS) by group with FN/FP at chosen cutoff
- Final table: `GTNA_predictions.txt` with `P_HS_supervised`, `P_HS_mixture`, etc.

## Reference

Code will be described in more detail upon submission of the associated paper.

