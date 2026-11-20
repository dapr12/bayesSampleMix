fit_small <- function(family = "gaussian", seed = 1) {
  d <- simulate_mixture_data(family, n_prob = 60, n_nonprob = 50, contamination = 0.2,
                             seed = seed)
  list(d = d, fit = mixture_fit(y ~ x1, d$prob, d$nonprob, family = family,
                                n_iter = 200, burnin = 50, thin = 2, seed = 1))
}

test_that("fitted object has the documented structure", {
  o <- fit_small()
  f <- o$fit
  expect_s3_class(f, "bsm_fit")
  expect_equal(nrow(f$draws), 75L)
  expect_equal(colnames(f$draws),
               c("(Intercept)", "x1", "sigma2", "contamination", "contamination_realized"))
  expect_equal(colnames(f$bad_draws), c("bad_mu0", "bad_tau2"))
  expect_length(f$filter_score, 50L)
  expect_true(all(f$filter_score >= 0 & f$filter_score <= 1))
  expect_true(all(f$indicator_entropy >= 0 & f$indicator_entropy <= log(2) + 1e-12))
  expect_equal(f$leverage, "x1")
  expect_equal(f$c_infl, 4)
  expect_true(all(f$draws[, "sigma2"] > 0))
  expect_true(all(f$draws[, "contamination"] > 0 & f$draws[, "contamination"] < 1))
})

test_that("coef, vcov, nobs, as.matrix, filter_scores and contamination agree with draws", {
  f <- fit_small()$fit
  expect_equal(coef(f), colMeans(f$draws[, 1:2]))
  expect_equal(coef(f, type = "median"), apply(f$draws[, 1:2], 2, median))
  expect_equal(vcov(f), cov(f$draws[, 1:2]))
  expect_equal(nobs(f), 110)
  expect_identical(as.matrix(f), f$draws)
  expect_equal(ncol(as.matrix(f, "all")), ncol(f$draws) + ncol(f$bad_draws))
  fs <- filter_scores(f, entropy = TRUE)
  expect_s3_class(fs, "data.frame")
  expect_equal(fs$filter_score, f$filter_score)
  ct <- contamination(f)
  expect_equal(ct[["mean"]], mean(f$draws[, "contamination"]))
  expect_true(ct[["lower"]] <= ct[["median"]] && ct[["median"]] <= ct[["upper"]])
})

test_that("summary and print methods produce output for all families", {
  for (fam in c("gaussian", "binomial", "poisson")) {
    f <- fit_small(fam)$fit
    s <- summary(f)
    expect_s3_class(s, "summary.bsm_fit")
    expect_equal(rownames(s$coefficients), c("(Intercept)", "x1"))
    expect_equal(s$coefficients[, "mean"], coef(f))
    expect_output(print(s), "Target-model coefficients")
    expect_output(print(f), "Contamination")
    if (fam == "poisson") expect_output(print(s), "acceptance")
    if (fam == "gaussian") expect_false(is.null(s$dispersion)) else expect_null(s$dispersion)
  }
})

test_that("predict() matches manual computation from the draws", {
  o <- fit_small("binomial")
  f <- o$fit
  nd <- data.frame(x1 = c(-1, 0, 2))
  p <- predict(f, nd)
  X <- cbind(1, nd$x1)
  manual <- plogis(f$draws[, 1:2] %*% t(X))
  expect_equal(p$fit, unname(colMeans(manual)))
  expect_equal(p$lwr, unname(apply(manual, 2, quantile, 0.025)))
  expect_true(all(p$lwr <= p$fit & p$fit <= p$upr))
  pl <- predict(f, nd, type = "link", interval = "none")
  expect_equal(names(pl), "fit")
  expect_equal(pl$fit, unname(colMeans(f$draws[, 1:2] %*% t(X))))
  dr <- predict(f, nd, summary = FALSE)
  expect_equal(dim(dr), c(nrow(f$draws), 3L))
  expect_error(predict(f), "newdata")
})

test_that("factor predictors use the union of levels and predict() handles them", {
  set.seed(2)
  dp <- data.frame(y = rnorm(60), x = rnorm(60), g = sample(c("a", "b"), 60, TRUE))
  dn <- data.frame(y = rnorm(40), x = rnorm(40), g = sample(c("a", "b", "c"), 40, TRUE))
  ## level "c" is absent from P, so the P design is rank deficient
  expect_error(mixture_fit(y ~ x + g, dp, dn, n_iter = 50, burnin = 10), "rank deficient")
  dn$g[dn$g == "c"] <- "a"
  f <- mixture_fit(y ~ x + g, dp, dn, n_iter = 100, burnin = 20, seed = 1)
  expect_equal(f$coef_names, c("(Intercept)", "x", "gb"))
  expect_equal(f$leverage, "x")
  p <- predict(f, data.frame(x = 0, g = c("a", "b")))
  expect_equal(p$fit[2] - p$fit[1], unname(coef(f)["gb"]))
})

test_that("plot() runs for every type without error", {
  f <- fit_small()$fit
  pdf(NULL)
  on.exit(dev.off())
  for (tp in c("trace", "density", "scores", "contamination")) {
    expect_invisible(plot(f, type = tp))
  }
  expect_invisible(plot(f, pars = c("sigma2", "bad_mu0")))
  expect_error(plot(f, pars = "nope"), "Unknown")
})

test_that("baseline_fit methods work", {
  d <- simulate_mixture_data("poisson", n_prob = 80, n_nonprob = 10, seed = 3)
  b <- baseline_fit(y ~ x1, d$prob, family = "poisson", n_iter = 200, burnin = 50, seed = 1)
  expect_s3_class(b, "bsm_baseline")
  expect_length(coef(b), 2L)
  expect_output(print(summary(b)), "Posterior summaries")
  expect_output(print(b), "Single-sample")
  expect_equal(nrow(predict(b, data.frame(x1 = 0:1))), 2L)
  expect_true(b$acceptance > 0 && b$acceptance <= 1)
})
