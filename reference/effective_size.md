# Effective sample size of MCMC draws

Computes the effective sample size of each column of a matrix of draws
from a single chain, using Geyer's (1992) initial monotone sequence
estimator of the integrated autocorrelation time.

## Usage

``` r
effective_size(x)
```

## Arguments

- x:

  Numeric vector or matrix (one column per parameter).

## Value

Named numeric vector of effective sample sizes. Constant columns return
`NA`.

## References

Geyer, C. J. (1992). Practical Markov chain Monte Carlo. *Statistical
Science*, 7(4), 473-483.

## Examples

``` r
set.seed(1)
x <- as.numeric(arima.sim(list(ar = 0.5), n = 2000))
effective_size(x)   # theory: about 2000 * (1 - 0.5) / (1 + 0.5) = 667
#> [1] 613.2262
```
