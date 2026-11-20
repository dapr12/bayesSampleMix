# bayesSampleMix 0.1.0

* Initial CRAN submission.

First packaged release of the contamination-mixture research code
(`mixture_filter.R`, `mixture_filter_glm.R`, `mixture_filter_pois.R`).

## Preserved from the original code

* Gaussian, logistic (Polya-Gamma) and Poisson (Laplace-proposal MH)
  samplers, models M0, Mjoint and M1, priors, initial values and default
  settings. Seeded fits reproduce the original functions draw for draw
  (verified by `tests/testthat/test-equivalence.R`).
* Baseline single-sample regressions (`baseline_fit()`), overlap diagnostics
  (`overlap_diagnostics()`, formerly `stage0()`), filter-score AUC
  (`filter_auc()`, formerly `auc_score()`).

## Changes

* M1: the rival regression now uses its exact full conditional whenever at
  least one unit is allocated to it (`mixture_control(m1_update = "exact")`).
  The original rule (prior draw when <= p + 1 units) is available with
  `m1_update = "original"`.
* All retained draws are stored, including sigma^2 and the contaminating
  component's parameters; realised contamination fraction added.
* Thinning, formula interface, factor handling, missing-data policy,
  input validation, seed handling that restores the global RNG state.
* `overlap_diagnostics()`/`filter_auc()` compute the AUC with the equivalent
  rank formula instead of an `n1 x n0` matrix.

## New

* S3 methods: `print()`, `summary()` (with ESS and MCSE), `coef()`, `vcov()`,
  `nobs()`, `predict()`, `plot()`, `as.matrix()`.
* `simulate_mixture_data()`, `effective_size()`, `filter_scores()`,
  `contamination()`.
