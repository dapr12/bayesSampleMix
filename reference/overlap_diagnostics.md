# Covariate-overlap diagnostics between the two samples

Pre-fitting ("Stage 0") diagnostics of how well the non-probability
sample overlaps the probability sample in the space of continuous
covariates. This is the `stage0()` function of the original code.

## Usage

``` r
overlap_diagnostics(prob, nonprob, vars = NULL)
```

## Arguments

- prob, nonprob:

  Numeric matrices or data frames of the continuous covariates for the P
  and NP samples (with matching column names).

- vars:

  Optional character vector of column names to use; by default all
  columns (which must be numeric).

## Value

A named list with elements `OC`, `EF` and `AUC_PS`, of class
`bsm_overlap`.

## Details

Mahalanobis distances of all units are computed with respect to the
probability-sample mean and covariance. Reported are

- `OC`:

  an overlap coefficient: the area shared by kernel density estimates of
  the P and NP distance distributions (on a 400-point grid), relative to
  the area of the P density. Values near 1 indicate good overlap.

- `EF`:

  the exceedance fraction: share of NP units whose distance exceeds the
  99th percentile of the P distances (about 0.01 under perfect overlap).

- `AUC_PS`:

  the area under the ROC curve of a logistic propensity model for sample
  membership (NP vs P) on the covariates; 0.5 indicates
  indistinguishable samples, values near 1 strong separation.

## Examples

``` r
set.seed(3)
P <- data.frame(x1 = rnorm(200), x2 = rnorm(200))
NP <- data.frame(x1 = rnorm(150, 0.5), x2 = rnorm(150))
overlap_diagnostics(P, NP)
#> Covariate overlap diagnostics (P vs NP)
#>   Overlap coefficient (OC):       0.942
#>   Exceedance fraction (EF, >q99): 0.020
#>   Propensity-score AUC:           0.656
```
