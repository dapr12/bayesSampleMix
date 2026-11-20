# Fit the contamination-mixture model from design matrices

Low-level interface that takes design matrices directly, mirroring the
original research functions `fit_mixture()`, `fit_logit_mixture()` and
`fit_pois_mixture()`. Most users should call
[`mixture_fit()`](https://dapr12.github.io/bayesSampleMix/reference/mixture_fit.md),
which builds the design matrices from a formula and calls this function.

## Usage

``` r
mixture_fit_matrix(
  Xp,
  yp,
  Xnp,
  ynp,
  family = "gaussian",
  model = c("Mjoint", "M0", "M1"),
  leverage = NULL,
  n_iter = 3000,
  burnin = 1500,
  thin = 1,
  seed = NULL,
  prior = mixture_prior(),
  control = mixture_control(),
  verbose = FALSE
)
```

## Arguments

- Xp, Xnp:

  Numeric design matrices (including an intercept column, if wanted) for
  the probability (P) and non-probability (NP) samples. Both must have
  the same columns.

- yp, ynp:

  Outcome vectors for the P and NP samples.

- family:

  One of `"gaussian"`, `"binomial"` (logistic regression for a 0/1,
  logical or two-level factor response) or `"poisson"`. A `family`
  object such as [`binomial()`](https://rdrr.io/r/stats/family.html) is
  also accepted (only its name is used; the links are fixed to identity,
  logit and log respectively).

- model:

  Contamination model: `"Mjoint"` (default), `"M0"` or `"M1"`; see the
  Model section.

- leverage:

  Integer column indices (or column names) of `Xp`/`Xnp` that enter the
  covariate density of model `"Mjoint"`; typically the continuous
  predictors. Equivalent to `cont_idx` in the original code. Ignored for
  models `"M0"` and `"M1"`.

- n_iter:

  Total number of MCMC iterations, including burn-in.

- burnin:

  Number of initial iterations discarded.

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

- verbose:

  If `TRUE`, report progress every 10% of iterations.

## Value

An object of class `bsm_fit`; see
[`mixture_fit()`](https://dapr12.github.io/bayesSampleMix/reference/mixture_fit.md).

## Examples

``` r
set.seed(1)
Xp <- cbind(1, rnorm(100)); yp <- drop(Xp %*% c(2, 5)) + rnorm(100)
Xnp <- cbind(1, rnorm(100)); ynp <- drop(Xnp %*% c(2, 5)) + rnorm(100)
ynp[1:25] <- rnorm(25, -20)          # 25% contaminated outcomes
fit <- mixture_fit_matrix(Xp, yp, Xnp, ynp, leverage = 2,
                          n_iter = 600, burnin = 200, seed = 1)
coef(fit)
#>       x1       x2 
#> 1.967023 5.059751 
```
