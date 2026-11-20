# Plot a contamination-mixture fit

Base-graphics displays of the posterior draws and filter scores.

## Usage

``` r
# S3 method for class 'bsm_fit'
plot(
  x,
  type = c("trace", "density", "scores", "contamination"),
  pars = NULL,
  ...
)
```

## Arguments

- x:

  A `bsm_fit` object.

- type:

  `"trace"` (trace plots of the retained draws), `"density"` (kernel
  density estimates of the marginal posteriors), `"scores"` (histogram
  of NP filter scores) or `"contamination"` (trace and histogram of
  \\1 - \pi\\).

- pars:

  Character vector of parameter names (columns of `as.matrix(x)`) for
  `"trace"` and `"density"`. Defaults to the target-model coefficients
  plus `contamination`.

- ...:

  Further graphical arguments passed to the underlying plotting
  functions.

## Value

`x`, invisibly.

## Examples

``` r
d <- simulate_mixture_data("gaussian", n_prob = 100, n_nonprob = 100,
                           contamination = 0.3, seed = 8)
fit <- mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 600, burnin = 200, seed = 1)
plot(fit, type = "trace")

plot(fit, type = "scores")
```
