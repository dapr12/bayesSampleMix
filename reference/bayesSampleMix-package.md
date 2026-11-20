# bayesSampleMix: Contamination-aware Bayesian integration of probability and non-probability samples

Fits Gaussian, logistic and Poisson regression models to a trusted
probability sample combined with a non-probability sample in which an
unknown subset of units may be incompatible with the target model. A
latent compatibility indicator per non-probability unit is sampled by
data augmentation, giving posterior inference for the regression
coefficients, the contamination fraction and unit-level filter scores.

## Details

Main functions:
[`mixture_fit()`](https://dapr12.github.io/bayesSampleMix/reference/mixture_fit.md),
[`baseline_fit()`](https://dapr12.github.io/bayesSampleMix/reference/baseline_fit.md),
[`overlap_diagnostics()`](https://dapr12.github.io/bayesSampleMix/reference/overlap_diagnostics.md),
[`simulate_mixture_data()`](https://dapr12.github.io/bayesSampleMix/reference/simulate_mixture_data.md).

## See also

Useful links:

- <https://github.com/dapr12/bayesSampleMix>

- Report bugs at <https://github.com/dapr12/bayesSampleMix/issues>

## Author

**Maintainer**: Diego Andres Perez Ruiz
<diego.perezruiz@manchester.ac.uk>

Authors:

- Constanza Avalos Valdebenito
  <constanza.avalos@postgrad.manchester.ac.uk>
