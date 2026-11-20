# Area under the ROC curve of filter scores against known labels

For simulation studies in which the true compatibility labels are known,
measures how well the filter scores rank compatible (label 1) above
contaminated (label 0) units. Ties count one half.

## Usage

``` r
filter_auc(labels, scores)
```

## Arguments

- labels:

  Vector of true labels (1 = compatible, 0 = contaminated).

- scores:

  Numeric scores, e.g.
  [`filter_scores()`](https://dapr12.github.io/bayesSampleMix/reference/filter_scores.md)
  of a fit.

## Value

A single number in \[0, 1\], or `NA` if either class is absent.

## Examples

``` r
filter_auc(c(1, 1, 0, 0), c(0.9, 0.7, 0.8, 0.1))
#> [1] 0.75
```
