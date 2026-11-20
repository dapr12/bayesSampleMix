# Single-sample Bayesian regression (reference estimators)

Fits the target regression model alone, with the same priors and the
same MCMC updates as the target component of
[`mixture_fit()`](https://dapr12.github.io/bayesSampleMix/reference/mixture_fit.md),
but without any contamination component. It reproduces the original
baseline functions `fit_bayes_lm()`, `fit_bayes_logit()` and
`fit_bayes_pois()` and is useful to compare

- the *probability-only* estimator (`data = prob_data`), and

- the *naively pooled* estimator
  (`data = rbind(prob_data, nonprob_data)`), which treats every NP unit
  as compatible,

with the contamination-aware fit.

## Usage

``` r
baseline_fit(
  formula,
  data,
  family = c("gaussian", "binomial", "poisson"),
  n_iter = 2000,
  burnin = 1000,
  thin = 1,
  seed = NULL,
  prior = mixture_prior(),
  control = mixture_control(),
  na_action = c("omit", "fail")
)

baseline_fit_matrix(
  x,
  y,
  family = "gaussian",
  n_iter = 2000,
  burnin = 1000,
  thin = 1,
  seed = NULL,
  prior = mixture_prior(),
  control = mixture_control()
)

# S3 method for class 'bsm_baseline'
summary(object, level = 0.95, ...)

# S3 method for class 'summary.bsm_baseline'
print(x, digits = 3, ...)

# S3 method for class 'bsm_baseline'
print(x, digits = 3, ...)

# S3 method for class 'bsm_baseline'
coef(object, ...)
```

## Arguments

- formula:

  A two-sided model formula, e.g. `y ~ x1 + x2`. Offsets are not
  supported.

- data:

  Data frame containing the variables in `formula`.

- family:

  One of `"gaussian"`, `"binomial"` (logistic regression for a 0/1,
  logical or two-level factor response) or `"poisson"`. A `family`
  object such as [`binomial()`](https://rdrr.io/r/stats/family.html) is
  also accepted (only its name is used; the links are fixed to identity,
  logit and log respectively).

- n_iter, burnin:

  Total iterations (including burn-in) and burn-in; the defaults
  (2000, 1000) are those of the original baseline functions.

- thin:

  Keep every `thin`-th post-burn-in draw.

- seed:

  Optional random seed. When supplied the fit is reproducible and the
  global random-number stream is restored afterwards; when `NULL` the
  current stream is used.

- prior:

  Prior hyperparameters, see
  [`mixture_prior()`](https://dapr12.github.io/bayesSampleMix/reference/mixture_prior.md).

- control:

  Algorithm settings, see
  [`mixture_control()`](https://dapr12.github.io/bayesSampleMix/reference/mixture_control.md).

- na_action:

  `"omit"` (default) drops rows with missing values in any model
  variable, separately within each sample, and records how many were
  dropped; `"fail"` stops with an error if any are present.

- x:

  For `baseline_fit_matrix()`, the design matrix; for the `print`
  methods, a `bsm_baseline` or `summary.bsm_baseline` object.

- y:

  Response vector (matrix interface).

- object:

  A `bsm_baseline` object.

- level:

  Credible-interval level for
  [`summary()`](https://rdrr.io/r/base/summary.html).

- ...:

  Not used.

- digits:

  Significant digits for printing.

## Value

An object of class `bsm_baseline` with elements `draws`, `coef_names`,
`family`, `acceptance` (Poisson only), `mcmc` and the information needed
by [`predict()`](https://rdrr.io/r/stats/predict.html).

## Examples

``` r
d <- simulate_mixture_data("gaussian", n_prob = 80, n_nonprob = 80,
                           contamination = 0.4, seed = 2)
prob_only <- baseline_fit(y ~ x1, d$prob, n_iter = 600, burnin = 200, seed = 1)
pooled <- baseline_fit(y ~ x1, rbind(d$prob, d$nonprob[names(d$prob)]),
                       n_iter = 600, burnin = 200, seed = 1)
rbind(prob_only = coef(prob_only), pooled = coef(pooled))
#>           (Intercept)       x1
#> prob_only    2.077655 4.885738
#> pooled       1.697441 2.774794
```
