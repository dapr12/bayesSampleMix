# Prior hyperparameters for the contamination-mixture model

Collects the hyperparameters of all prior distributions used by
[`mixture_fit()`](https://dapr12.github.io/bayesSampleMix/reference/mixture_fit.md)
and
[`baseline_fit()`](https://dapr12.github.io/bayesSampleMix/reference/baseline_fit.md).
The defaults are exactly those of the original research code, so that
results are reproducible.

## Usage

``` r
mixture_prior(
  beta_scale = 25,
  sigma_shape = 2,
  sigma_rate = 1,
  pi_shape1 = 1,
  pi_shape2 = 1,
  tau_shape = 2,
  tau_rate = 2,
  mu0_var = 100,
  q_shape1 = 1,
  q_shape2 = 1,
  lambda_shape = 1,
  lambda_rate = 1
)
```

## Arguments

- beta_scale:

  Prior variance \\b_0\\ of each regression coefficient: \\\beta \sim
  N_p(0, b_0 I_p)\\. Also used for the rival regression \\\beta\_{bad}\\
  of model `"M1"`.

- sigma_shape, sigma_rate:

  Inverse-gamma prior for the Gaussian error variance, \\\sigma^2 \sim
  IG(a\_\sigma, b\_\sigma)\\ (Gaussian family only; also used for the
  rival-regression variance in Gaussian `"M1"`).

- pi_shape1, pi_shape2:

  Beta prior for the probability \\\pi\\ that a non-probability
  observation is compatible ("good"): \\\pi \sim Beta(a\_\pi, b\_\pi)\\.
  The contamination fraction is \\1 - \pi\\.

- tau_shape, tau_rate:

  Inverse-gamma prior for the variance \\\tau^2\\ of the covariate-free
  contaminating outcome density (Gaussian `"M0"` and `"Mjoint"`).

- mu0_var:

  Prior variance of the mean \\\mu_0\\ of the contaminating outcome
  density (Gaussian `"M0"`/`"Mjoint"`). The prior mean is the sample
  mean of the non-probability outcomes (an empirical-Bayes choice
  inherited from the original code).

- q_shape1, q_shape2:

  Beta prior for the base rate \\q\\ of the contaminating Bernoulli
  density (binomial `"M0"`/`"Mjoint"`).

- lambda_shape, lambda_rate:

  Gamma prior (shape, rate) for the rate \\\lambda_0\\ of the
  contaminating Poisson density (Poisson `"M0"`/`"Mjoint"`).

## Value

An object of class `bsm_prior` (a named list).

## See also

[`mixture_fit()`](https://dapr12.github.io/bayesSampleMix/reference/mixture_fit.md),
[`mixture_control()`](https://dapr12.github.io/bayesSampleMix/reference/mixture_control.md)

## Examples

``` r
mixture_prior()
#> Prior hyperparameters (bayesSampleMix)
#>   beta_scale  sigma_shape   sigma_rate    pi_shape1    pi_shape2    tau_shape 
#>           25            2            1            1            1            2 
#>     tau_rate      mu0_var     q_shape1     q_shape2 lambda_shape  lambda_rate 
#>            2          100            1            1            1            1 
mixture_prior(pi_shape1 = 2, pi_shape2 = 1)  # a priori favour compatibility
#> Prior hyperparameters (bayesSampleMix)
#>   beta_scale  sigma_shape   sigma_rate    pi_shape1    pi_shape2    tau_shape 
#>           25            2            1            2            1            2 
#>     tau_rate      mu0_var     q_shape1     q_shape2 lambda_shape  lambda_rate 
#>            2          100            1            1            1            1 
```
