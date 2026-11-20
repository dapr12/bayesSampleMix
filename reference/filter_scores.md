# Filter scores and contamination summaries

`filter_scores()` returns the posterior probability \\s_j\\ that each
non-probability unit is compatible with the target model (the
Rao-Blackwellised posterior mean of its latent indicator). Small values
flag likely contaminated units. `contamination()` summarises the
posterior of the contamination fraction \\1 - \pi\\.

## Usage

``` r
filter_scores(object, entropy = FALSE)

contamination(object, level = 0.95)
```

## Arguments

- object:

  A `bsm_fit` object.

- entropy:

  If `TRUE`, return a data frame that also contains the binary entropy
  (nats) of each score, a measure of classification uncertainty
  (maximal, \\\log 2\\, at \\s_j = 0.5\\).

- level:

  Credible-interval level.

## Value

`filter_scores()`: a numeric vector named by the rows of `nonprob_data`
used in the fit, or a data frame if `entropy = TRUE`. `contamination()`:
a named vector with posterior mean, median and credible limits of \\1 -
\pi\\.

## Examples

``` r
d <- simulate_mixture_data("gaussian", n_prob = 100, n_nonprob = 100,
                           contamination = 0.3, seed = 6)
fit <- mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 600, burnin = 200, seed = 1)
s <- filter_scores(fit)
tapply(s, d$nonprob$contaminated, mean)   # low scores for contaminated units
#>         0         1 
#> 0.9486064 0.2355229 
contamination(fit)
#>      mean    median     lower     upper 
#> 0.2690224 0.2678621 0.1862934 0.3651758 
```
