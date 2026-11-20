# Posterior summaries of a contamination-mixture fit

Posterior summaries of a contamination-mixture fit

## Usage

``` r
# S3 method for class 'bsm_fit'
summary(object, level = 0.95, ...)

# S3 method for class 'summary.bsm_fit'
print(x, digits = 3, ...)

# S3 method for class 'bsm_fit'
print(x, digits = 3, ...)

# S3 method for class 'bsm_fit'
coef(object, type = c("mean", "median"), ...)

# S3 method for class 'bsm_fit'
vcov(object, ...)

# S3 method for class 'bsm_fit'
nobs(object, ...)
```

## Arguments

- object, x:

  A `bsm_fit` object from
  [`mixture_fit()`](https://dapr12.github.io/bayesSampleMix/reference/mixture_fit.md)
  or
  [`mixture_fit_matrix()`](https://dapr12.github.io/bayesSampleMix/reference/mixture_fit_matrix.md)
  (or a `bsm_baseline` object for the methods that also apply to
  [`baseline_fit()`](https://dapr12.github.io/bayesSampleMix/reference/baseline_fit.md)).

- level:

  Probability mass of the equal-tailed credible intervals.

- ...:

  Not used.

- digits:

  Number of significant digits to print.

- type:

  For [`coef()`](https://rdrr.io/r/stats/coef.html): posterior `"mean"`
  (default) or `"median"`.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.bsm_fit` with elements

- coefficients:

  posterior mean, standard deviation, credible interval, median,
  effective sample size and Monte Carlo standard error of each
  target-model coefficient;

- dispersion:

  the same for \\\sigma^2\\ (Gaussian family only);

- contamination:

  the same for \\1 - \pi\\ and the realised contamination fraction;

- filter:

  number of NP units, posterior expected number of compatible units
  \\\sum_j s_j\\, number of units with \\s_j \< 0.5\\, and mean
  indicator entropy;

plus model settings. All quantities are computed from the retained
posterior draws; no asymptotic approximations are used.

## Examples

``` r
d <- simulate_mixture_data("gaussian", n_prob = 100, n_nonprob = 100,
                           contamination = 0.25, seed = 4)
fit <- mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 600, burnin = 200, seed = 1)
summary(fit)
#> Bayesian contamination-mixture fit (bayesSampleMix)
#> Call: mixture_fit(formula = y ~ x1, prob_data = d$prob, nonprob_data = d$nonprob, 
#>     n_iter = 600, burnin = 200, seed = 1)
#> Family:  Gaussian (identity link)
#> Model:   Mjoint (contaminating density + covariate leverage term)
#> Leverage covariates: x1  (inflation c = 4)
#> Samples: n_P = 100, n_NP = 100
#> MCMC:    600 iterations, burn-in 200, thin 1 -> 400 draws
#> 
#> Target-model coefficients (posterior; 95% credible intervals):
#>             mean     sd 2.5%  50% 97.5% ess    mcse
#> (Intercept) 1.94 0.0786 1.78 1.94  2.10 349 0.00421
#> x1          4.91 0.0835 4.74 4.91  5.07 269 0.00509
#> 
#> Error variance:
#>         mean    sd  2.5%   50% 97.5% ess    mcse
#> sigma2 0.946 0.107 0.747 0.945  1.17 322 0.00595
#> 
#> Contamination:
#>                         mean     sd  2.5%   50% 97.5% ess    mcse
#> contamination          0.245 0.0506 0.159 0.246 0.351 216 0.00344
#> contamination_realized 0.241 0.0248 0.200 0.240 0.290 182 0.00184
#> 
#> Filter scores of NP units (s_j = posterior P(compatible)):
#>   expected compatible units: 76.0 of 100; units with s_j < 0.5: 20; mean entropy: 0.167
coef(fit)
#> (Intercept)          x1 
#>    1.943749    4.911094 
contamination(fit)
#>      mean    median     lower     upper 
#> 0.2452251 0.2457701 0.1587920 0.3513224 
```
