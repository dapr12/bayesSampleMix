# Posterior predictions from the target regression model

Computes, for each row of `newdata`, the posterior distribution of the
target-model linear predictor (`type = "link"`) or mean response
(`type = "response"`: \\x'\beta\\, \\\mathrm{logit}^{-1}(x'\beta)\\ or
\\\exp(x'\beta)\\), using every retained posterior draw. This is the
regression function of the population represented by the probability
sample and the compatible NP units; it does not include
observation-level noise.

## Usage

``` r
# S3 method for class 'bsm_fit'
predict(
  object,
  newdata,
  type = c("response", "link"),
  interval = c("credible", "none"),
  level = 0.95,
  summary = TRUE,
  ...
)

# S3 method for class 'bsm_baseline'
predict(
  object,
  newdata,
  type = c("response", "link"),
  interval = c("credible", "none"),
  level = 0.95,
  summary = TRUE,
  ...
)
```

## Arguments

- object:

  A `bsm_fit` or `bsm_baseline` object.

- newdata:

  Data frame (formula fits) or design matrix (matrix fits) at which to
  predict. Required.

- type:

  `"response"` (default) or `"link"`.

- interval:

  `"credible"` (default) to add equal-tailed credible limits, or
  `"none"`.

- level:

  Credible-interval level.

- summary:

  If `FALSE`, return the full matrix of posterior draws (draws in rows,
  observations in columns) instead of summaries.

- ...:

  Not used.

## Value

A data frame with columns `fit` (posterior mean) and, if requested,
`lwr` and `upr`; or a matrix of draws when `summary = FALSE`.

## Examples

``` r
d <- simulate_mixture_data("binomial", n_prob = 150, n_nonprob = 150,
                           contamination = 0.3, seed = 7)
fit <- mixture_fit(y ~ x1, d$prob, d$nonprob, family = "binomial",
                   n_iter = 500, burnin = 200, seed = 1)
predict(fit, newdata = data.frame(x1 = c(-1, 0, 1)))
#>         fit       lwr       upr
#> 1 0.3760308 0.2326381 0.5099556
#> 2 0.5514263 0.4528796 0.6349134
#> 3 0.7154087 0.6233067 0.8109616
```
