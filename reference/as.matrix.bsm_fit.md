# Posterior draws as a matrix

Posterior draws as a matrix

## Usage

``` r
# S3 method for class 'bsm_fit'
as.matrix(x, which = c("target", "bad", "all"), ...)

# S3 method for class 'bsm_baseline'
as.matrix(x, ...)
```

## Arguments

- x:

  A `bsm_fit` or `bsm_baseline` object.

- which:

  For `bsm_fit`: `"target"` (target-model parameters and contamination
  quantities, default), `"bad"` (contaminating-component parameters) or
  `"all"`.

- ...:

  Not used.

## Value

A numeric matrix with one row per retained draw; suitable for
[`coda::mcmc()`](https://rdrr.io/pkg/coda/man/mcmc.html).

## Examples

``` r
d <- simulate_mixture_data("poisson", n_prob = 80, n_nonprob = 80,
                           contamination = 0.2, seed = 5)
fit <- mixture_fit(y ~ x1, d$prob, d$nonprob, family = "poisson",
                   n_iter = 400, burnin = 100, seed = 1)
head(as.matrix(fit))
#>      (Intercept)        x1 contamination contamination_realized
#> [1,]    1.221602 0.4297103     0.3450401                 0.3250
#> [2,]    1.332617 0.3822226     0.2808199                 0.2750
#> [3,]    1.261872 0.3790944     0.3287653                 0.3125
#> [4,]    1.303046 0.3674293     0.2247163                 0.3125
#> [5,]    1.322738 0.4710368     0.2858494                 0.2375
#> [6,]    1.245842 0.4930108     0.2849340                 0.3125
```
