# Algorithmic settings for the contamination-mixture samplers

Algorithmic settings for the contamination-mixture samplers

## Usage

``` r
mixture_control(
  c_infl = NULL,
  pg_terms = 128L,
  irls_steps = 3L,
  m1_update = c("exact", "original")
)
```

## Arguments

- c_infl:

  Variance-inflation factor \\c \> 1\\ of the covariate density of
  contaminated observations in model `"Mjoint"`: good observations have
  \\x \sim N(\hat\mu_P, \hat\Sigma_P)\\, contaminated ones \\x \sim
  N(\hat\mu_P, c\\\hat\Sigma_P)\\. `NULL` (default) uses the
  family-specific default of the original code: 4 for Gaussian, 2.5 for
  binomial and Poisson.

- pg_terms:

  Number of terms \\K\\ in the truncated sum-of-gammas representation
  used to draw Polya-Gamma variables (binomial family). Default 128, as
  in the original code.

- irls_steps:

  Number of Newton (IRLS) steps used to locate the conditional mode for
  the Laplace proposal of the Poisson coefficient update. Default 3, as
  in the original code. Larger values make the independence
  Metropolis-Hastings step closer to exact (see Details in
  [`mixture_fit()`](https://dapr12.github.io/bayesSampleMix/reference/mixture_fit.md)).

- m1_update:

  How the rival-regression component of model `"M1"` is updated when few
  observations are allocated to it. `"exact"` (default) always uses the
  correct full conditional given the allocated observations (a proper
  prior makes this well defined for any number of observations).
  `"original"` reproduces the original code, which draws
  \\\beta\_{bad}\\ from its prior whenever \\\le p + 1\\ observations
  are allocated to the contaminated component; this ignores those
  observations and is not a valid Gibbs step, so it is provided only for
  reproducing earlier results.

## Value

An object of class `bsm_control`.

## Examples

``` r
mixture_control(c_infl = 3)
#> $c_infl
#> [1] 3
#> 
#> $pg_terms
#> [1] 128
#> 
#> $irls_steps
#> [1] 3
#> 
#> $m1_update
#> [1] "exact"
#> 
#> attr(,"class")
#> [1] "bsm_control"
```
