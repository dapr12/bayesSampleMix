# Contamination-aware Bayesian integration of probability and non-probability samples

Fits a regression model to a trusted probability (P) sample and a
non-probability (NP) sample in which an unknown subset of observations
may be incompatible with the target model ("contaminated"). Each NP
observation carries a latent indicator of compatibility; compatible
observations are pooled with the P sample to estimate the target
regression, while contaminated observations are absorbed by a separate
component.

## Usage

``` r
mixture_fit(
  formula,
  prob_data,
  nonprob_data,
  family = c("gaussian", "binomial", "poisson"),
  model = c("Mjoint", "M0", "M1"),
  leverage = NULL,
  n_iter = 3000,
  burnin = 1500,
  thin = 1,
  seed = NULL,
  prior = mixture_prior(),
  control = mixture_control(),
  na_action = c("omit", "fail"),
  verbose = FALSE
)
```

## Arguments

- formula:

  A two-sided model formula, e.g. `y ~ x1 + x2`. Offsets are not
  supported.

- prob_data:

  Data frame with the probability (trusted) sample.

- nonprob_data:

  Data frame with the non-probability sample. Must contain all variables
  used in `formula`.

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

  Names of the design-matrix columns entering the covariate density of
  `"Mjoint"`. `NULL` (default) selects all non-constant columns taking
  more than two distinct values in the P sample (i.e. continuous
  predictors; intercept and dummy variables are excluded).

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

- na_action:

  `"omit"` (default) drops rows with missing values in any model
  variable, separately within each sample, and records how many were
  dropped; `"fail"` stops with an error if any are present.

- verbose:

  If `TRUE`, report progress every 10% of iterations.

## Value

An object of class `bsm_fit`, a list including

- draws:

  matrix of retained posterior draws of the target-model coefficients
  (plus `sigma2` for Gaussian), `contamination` and
  `contamination_realized`;

- bad_draws:

  draws of the contaminating-component parameters;

- filter_score:

  posterior compatibility probabilities \\s_j\\ of the NP observations
  used in the fit;

- indicator_entropy:

  binary entropy (nats) of each \\s_j\\;

- nonprob_rows:

  row names of `nonprob_data` used in the fit (after removing missing
  values), aligned with `filter_score`;

and settings needed by the methods
[`summary.bsm_fit()`](https://dapr12.github.io/bayesSampleMix/reference/summary.bsm_fit.md),
[`predict.bsm_fit()`](https://dapr12.github.io/bayesSampleMix/reference/predict.bsm_fit.md),
[`plot.bsm_fit()`](https://dapr12.github.io/bayesSampleMix/reference/plot.bsm_fit.md),
[`coef.bsm_fit()`](https://dapr12.github.io/bayesSampleMix/reference/summary.bsm_fit.md)
and
[`as.matrix.bsm_fit()`](https://dapr12.github.io/bayesSampleMix/reference/as.matrix.bsm_fit.md).

## Model

For P-sample units \\i\\ the outcome follows the target model \\f(y_i
\mid x_i, \theta)\\ (Gaussian, logistic or Poisson regression with
coefficients \\\beta\\). For NP units \\j\\, \$\$z_j \sim
\mathrm{Bernoulli}(\pi),\$\$ \$\$p(y_j, x_j \mid z_j = 1) = f(y_j \mid
x_j, \theta)\\ g_1(x_j),\qquad p(y_j, x_j \mid z_j = 0) = h(y_j \mid
x_j, \phi)\\ g_0(x_j).\$\$ The contaminating outcome density \\h\\
depends on `model`:

- `"M0"`:

  a covariate-free density: \\N(\mu_0, \tau^2)\\ (Gaussian),
  \\Bernoulli(q)\\ (binomial), \\Poisson(\lambda_0)\\ (Poisson); and
  \\g_1 = g_0 = 1\\.

- `"Mjoint"`:

  as `"M0"`, plus a covariate ("leverage") term for the continuous
  predictors \\x_c\\ selected by `leverage`: \\g_1 = N(x_c; \hat\mu_P,
  \hat\Sigma_P)\\ and \\g_0 = N(x_c; \hat\mu_P, c\\\hat\Sigma_P)\\,
  where \\\hat\mu_P, \hat\Sigma_P\\ are the P-sample mean and covariance
  (fixed plug-in values) and \\c\\ is `control$c_infl`. NP points far
  from the P covariate distribution (in Mahalanobis distance) are
  thereby more likely to be flagged.

- `"M1"`:

  a rival regression of the same family with its own coefficients
  \\\beta\_{bad}\\ (and variance, Gaussian case); no covariate term.

Because the P sample contributes only to the target component, the
labels of the two mixture components are identified (no label
switching).

## Algorithm

A data-augmentation Gibbs sampler alternates between (1) drawing every
\\z_j\\ from its Bernoulli full conditional, (2) drawing \\\pi \sim
Beta(a\_\pi + \sum z_j, b\_\pi + n\_{NP} - \sum z_j)\\, (3) updating the
target-model parameters given the P sample plus the NP units with \\z_j
= 1\\, and (4) updating the contaminating component given the units with
\\z_j = 0\\. Step (3) is family specific and follows the original
implementation:

- Gaussian: conjugate normal / inverse-gamma Gibbs updates.

- Binomial: Polya-Gamma augmentation (Polson, Scott and Windle, 2013),
  with PG(1, z) variables drawn from a `pg_terms`-term truncation of
  their sum-of-gammas representation plus a tail-mean correction. This
  is a close approximation rather than an exact PG draw.

- Poisson: a Metropolis-Hastings step with a Laplace (Gaussian) proposal
  centred at the conditional mode, which is found by `irls_steps` Newton
  iterations from the current value. Once these iterations have
  converged the proposal is state independent and the step is an exact
  independence sampler; otherwise it is an approximation. The
  post-burn-in acceptance rate is reported by
  [`summary.bsm_fit()`](https://dapr12.github.io/bayesSampleMix/reference/summary.bsm_fit.md).

The filter score of NP unit \\j\\ is the Rao-Blackwellised posterior
probability of compatibility, \\s_j = E\[P(z_j = 1 \mid \cdot) \mid
data\]\\, averaged over all post-burn-in iterations (thinning affects
only the stored parameter draws).

## Contamination summaries

Two posterior quantities are reported: `contamination` is \\1 - \pi\\,
the population-level probability that an NP unit is contaminated (as in
the original code); `contamination_realized` is the fraction of the
observed NP units allocated to the contaminating component at each
iteration, \\1 - \sum_j z_j / n\_{NP}\\.

\\1 - \pi\\ is the weight of the contaminating component, which is
identified only to the extent that the two components can be told apart.
With the covariate-free contaminating densities of `"M0"`/`"Mjoint"`,
compatible units that fit both components about equally well are partly
allocated to the contaminating one, and the `"Mjoint"` covariate term
also allocates compatible units in the tails of the P covariate
distribution to it. In simulations with binomial and Poisson outcomes
this led to posterior means of \\1 - \pi\\ well above the true share of
contaminated units, while the target coefficients remained accurate (see
the package vignettes). Read \\1 - \pi\\ as a conservative measure of
incompatibility rather than a count of contaminated units.

## References

Polson, N. G., Scott, J. G. and Windle, J. (2013). Bayesian inference
for logistic models using Polya-Gamma latent variables. *Journal of the
American Statistical Association*, 108(504), 1339-1349.

Diebolt, J. and Robert, C. P. (1994). Estimation of finite mixture
distributions through Bayesian sampling. *Journal of the Royal
Statistical Society: Series B*, 56(2), 363-375.

## See also

[`baseline_fit()`](https://dapr12.github.io/bayesSampleMix/reference/baseline_fit.md)
for probability-only and naively pooled reference fits,
[`overlap_diagnostics()`](https://dapr12.github.io/bayesSampleMix/reference/overlap_diagnostics.md)
for covariate-overlap checks,
[`simulate_mixture_data()`](https://dapr12.github.io/bayesSampleMix/reference/simulate_mixture_data.md)
for example data.

## Examples

``` r
d <- simulate_mixture_data("gaussian", n_prob = 100, n_nonprob = 100,
                           contamination = 0.3, seed = 1)
fit <- mixture_fit(y ~ x1, d$prob, d$nonprob, family = "gaussian",
                   n_iter = 600, burnin = 200, seed = 42)
fit
#> Bayesian contamination-mixture fit (bayesSampleMix)
#> Family: gaussian   Model: Mjoint
#> n_P = 100, n_NP = 100, 400 posterior draws
#> 
#> Posterior mean coefficients:
#> (Intercept)          x1 
#>        2.08        4.97 
#> 
#> Contamination (1 - pi): 0.323  [95% CrI 0.224, 0.419]
summary(fit)
#> Bayesian contamination-mixture fit (bayesSampleMix)
#> Call: mixture_fit(formula = y ~ x1, prob_data = d$prob, nonprob_data = d$nonprob, 
#>     family = "gaussian", n_iter = 600, burnin = 200, seed = 42)
#> Family:  Gaussian (identity link)
#> Model:   Mjoint (contaminating density + covariate leverage term)
#> Leverage covariates: x1  (inflation c = 4)
#> Samples: n_P = 100, n_NP = 100
#> MCMC:    600 iterations, burn-in 200, thin 1 -> 400 draws
#> 
#> Target-model coefficients (posterior; 95% credible intervals):
#>             mean     sd 2.5%  50% 97.5% ess    mcse
#> (Intercept) 2.08 0.0792 1.92 2.08  2.22 357 0.00419
#> x1          4.97 0.0821 4.83 4.97  5.14 393 0.00414
#> 
#> Error variance:
#>         mean    sd  2.5%   50% 97.5% ess    mcse
#> sigma2 0.956 0.105 0.774 0.948  1.18 297 0.00608
#> 
#> Contamination:
#>                         mean     sd  2.5%   50% 97.5% ess    mcse
#> contamination          0.323 0.0515 0.224 0.321 0.419 230 0.00339
#> contamination_realized 0.316 0.0288 0.260 0.320 0.370 231 0.00189
#> 
#> Filter scores of NP units (s_j = posterior P(compatible)):
#>   expected compatible units: 68.6 of 100; units with s_j < 0.5: 26; mean entropy: 0.212
head(filter_scores(fit))
#>            1            2            3            4            5            6 
#> 3.329450e-04 5.697014e-35 9.349462e-01 9.238125e-01 4.716416e-17 4.039575e-06 
```
