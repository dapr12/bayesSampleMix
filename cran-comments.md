## Submission

This is the first submission of bayesSampleMix to CRAN.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Diego Andres Perez Ruiz <diego.perezruiz@manchester.ac.uk>'
  New submission

## Test environments

* local: R 4.6.2
* GitHub Actions: macOS-latest (release), windows-latest (release),
  ubuntu-latest (devel, release, oldrel-1)
* win-builder: R-devel

## Notes for reviewers

* The package depends only on base R (`stats`, `graphics`).
* The longer statistical-recovery tests are skipped on CRAN with
  `testthat::skip_on_cran()` to keep check time short. They run on
  GitHub Actions on every push.
* Functions that accept a `seed` argument restore the user's global RNG
  state on exit, so the global environment is not modified.
* Possibly misspelled words in DESCRIPTION (Polya, Gamma, Mahalanobis,
  Metropolis, Hastings) are names of statistical methods and are correct.
