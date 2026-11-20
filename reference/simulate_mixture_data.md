# Simulate probability and contaminated non-probability samples

Generates example data with a known data-generating process for
demonstrations, tests and simulation studies.

## Usage

``` r
simulate_mixture_data(
  family = c("gaussian", "binomial", "poisson"),
  n_prob = 200,
  n_nonprob = 200,
  contamination = 0.25,
  beta = NULL,
  sigma = 1,
  bad_beta = NULL,
  x_shift = 0,
  outcome_shift = 0,
  seed = NULL
)
```

## Arguments

- family:

  `"gaussian"`, `"binomial"` or `"poisson"`.

- n_prob, n_nonprob:

  Sample sizes.

- contamination:

  Fraction of contaminated NP units, in \[0, 1\].

- beta:

  True target coefficients (intercept first). Defaults: `c(2, 5)`
  (Gaussian), `c(0.2, 1, -0.8)` (binomial), `c(1.2, 0.5, -0.4)`
  (Poisson).

- sigma:

  Error standard deviation (Gaussian only).

- bad_beta:

  Coefficients of the contaminating mechanism. Default: the intercept of
  `beta` with all slopes sign-flipped.

- x_shift:

  Mean shift of contaminated units' covariates.

- outcome_shift:

  Constant added to contaminated units' linear predictor.

- seed:

  Optional seed (the global RNG state is restored afterwards).

## Value

A list with data frames `prob` (columns `y`, `x1`, ...) and `nonprob`
(same columns plus `contaminated`, 1 for contaminated units), and the
list `truth` of generating values.

## Details

Covariates \\x_1, \dots, x\_{k}\\ (\\k\\ = `length(beta) - 1`) are
independent standard normal in the probability sample and in the
compatible part of the non-probability sample. Exactly
`round(contamination * n_nonprob)` NP units are contaminated: their
covariates are drawn from \\N(x\_{shift}, 1)\\ and their outcomes from
the same family with linear predictor \\x'\beta\_{bad} + \delta\\
(`bad_beta`, `outcome_shift`). NP rows are randomly permuted.

## Examples

``` r
d <- simulate_mixture_data("poisson", n_prob = 50, n_nonprob = 40,
                           contamination = 0.25, seed = 1)
str(d$nonprob)
#> 'data.frame':    40 obs. of  4 variables:
#>  $ y           : int  1 1 2 1 1 1 1 4 2 1 ...
#>  $ x1          : num  0.0174 1.5866 -0.6817 -1.1159 -0.3309 ...
#>  $ x2          : num  2.075 -0.4 -1.384 -0.334 -1.37 ...
#>  $ contaminated: int  0 1 0 0 1 0 0 0 1 0 ...
table(d$nonprob$contaminated)
#> 
#>  0  1 
#> 30 10 
```
