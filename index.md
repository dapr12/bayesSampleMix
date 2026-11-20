# bayesSampleMix

**Contamination-aware Bayesian integration of probability and
non-probability samples.**

`bayesSampleMix` fits Gaussian, logistic and Poisson regressions to a
trusted probability sample (P) combined with a non-probability sample
(NP) in which an unknown subset of units may be incompatible with the
target model. Each NP unit has a latent compatibility indicator;
compatible units are pooled with the P sample, contaminated units are
absorbed by a separate mixture component. You get

- posterior inference for the target regression coefficients,
- the posterior of the contamination fraction,
- a **filter score** for every NP unit — its posterior probability of
  being compatible,
- covariate-overlap diagnostics and reference fits (probability-only and
  naively pooled).

## Why not just pool?

| Strategy | Bias under contamination | Precision |
|----|----|----|
| Probability sample only | none | low |
| Naive pooling (treat all NP units as valid) | yes — does not vanish as n grows | high (but around the wrong value) |
| **Contamination mixture** | adapts: uses only units the data deem compatible | between the two |

## Installation

Once the package is on CRAN:

``` r

install.packages("bayesSampleMix")
```

The development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("dapr12/bayesSampleMix", build_vignettes = TRUE)
```

The package needs only base R (`stats`, `graphics`); `knitr`,
`rmarkdown`, `testthat` and `coda` are suggested.

## Quick start

``` r

library(bayesSampleMix)

# simulated P and NP samples; 30% of NP units follow a different mechanism
d <- simulate_mixture_data("gaussian", n_prob = 150, n_nonprob = 300,
                           contamination = 0.3, beta = c(2, 5),
                           x_shift = 1.5, seed = 2026)

fit <- mixture_fit(
  formula      = y ~ x1,
  prob_data    = d$prob,
  nonprob_data = d$nonprob,
  family       = "gaussian",   # or "binomial", "poisson"
  model        = "Mjoint",     # or "M0", "M1"
  n_iter       = 3000,
  burnin       = 1000,
  thin         = 1,
  seed         = 123
)

summary(fit)                 # coefficients, contamination, filter scores, ESS
coef(fit)                    # posterior means
contamination(fit)           # posterior of 1 - pi
s <- filter_scores(fit)      # posterior P(compatible) for each NP unit
predict(fit, newdata = data.frame(x1 = c(-1, 0, 1)))
plot(fit, type = "trace")

# reference estimators
baseline_fit(y ~ x1, d$prob)                                   # probability only
baseline_fit(y ~ x1, rbind(d$prob, d$nonprob[names(d$prob)]))  # naive pooling
```

## Interpreting the results

- **Coefficients** are posterior means, standard deviations and
  equal-tailed credible intervals computed from the MCMC draws, with
  effective sample sizes and Monte Carlo standard errors. They describe
  the population represented by the P sample.
- **Contamination** `1 - pi` is the posterior probability that an NP
  unit is incompatible; `contamination_realized` is the share of the
  observed NP units allocated to the contaminating component.
- **Filter scores** close to 0 flag units that the data consider
  incompatible with the target model; scores near 0.5 (high indicator
  entropy) indicate units that cannot be classified confidently.
- Check
  [`overlap_diagnostics()`](https://dapr12.github.io/bayesSampleMix/reference/overlap_diagnostics.md)
  before fitting and the trace plots after.

## Models and algorithms

| `family` | target model | contaminating density (M0/Mjoint) | target update |
|----|----|----|----|
| `gaussian` | N(x’β, σ²) | N(μ0, τ²) | conjugate Gibbs |
| `binomial` | Bernoulli(logit⁻¹(x’β)) | Bernoulli(q) | Pólya–Gamma augmentation |
| `poisson` | Poisson(exp(x’β)) | Poisson(λ0) | Laplace-proposal Metropolis–Hastings |

`Mjoint` adds a covariate term that compares each NP unit’s continuous
covariates with the P-sample distribution (Mahalanobis-type); `M1`
replaces the contaminating density by a rival regression. See
[`?mixture_fit`](https://dapr12.github.io/bayesSampleMix/reference/mixture_fit.md)
and the vignettes:

``` r

vignette("introduction", package = "bayesSampleMix")
vignette("gaussian-model", package = "bayesSampleMix")
vignette("logistic-model", package = "bayesSampleMix")
vignette("poisson-model", package = "bayesSampleMix")
```

## Relation to the original research code

The samplers, priors, initial values and defaults of the original
scripts (`mixture_filter.R`, `mixture_filter_glm.R`,
`mixture_filter_pois.R`) are preserved; the test suite checks that
seeded fits reproduce them draw for draw. Differences and known
approximations are documented in `NEWS.md` (notably the corrected M1
update, available in its original form via
`mixture_control(m1_update = "original")`). The original Monte Carlo
studies are in `inst/validation/` and the rent-index application in
`inst/scripts/`.

## Authors

- **Diego Andrés Pérez Ruiz** (maintainer), Department of Social
  Statistics, University of Manchester
- **Constanza Avalos Valdebenito**, Department of Social Statistics,
  University of Manchester

## Citation

If you use `bayesSampleMix`, please cite the package and the methodology
paper. Run `citation("bayesSampleMix")` in R for BibTeX entries.

> Pérez Ruiz, D. A. and Avalos Valdebenito, C. (2026). *bayesSampleMix:
> Contamination-Aware Bayesian Integration of Probability and
> Non-Probability Samples*. R package version 0.1.0.
> <https://github.com/dapr12/bayesSampleMix>

> Pérez Ruiz, D. A. and Avalos Valdebenito, C. (2026). *A Fully Bayesian
> Contamination-Mixture Framework for Probability / Non-Probability Data
> Integration*. Manuscript.

## Contributing

Bug reports and pull requests are welcome at
<https://github.com/dapr12/bayesSampleMix/issues>. See the [contributing
guide](https://dapr12.github.io/bayesSampleMix/CONTRIBUTING.md). This
project follows a [Code of
Conduct](https://dapr12.github.io/bayesSampleMix/CODE_OF_CONDUCT.md).

## License

MIT © 2026 Diego Andrés Pérez Ruiz and Constanza Avalos Valdebenito. See
[LICENSE.md](https://dapr12.github.io/bayesSampleMix/LICENSE.md).
