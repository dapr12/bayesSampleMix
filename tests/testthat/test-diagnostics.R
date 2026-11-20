test_that("effective_size() is close to theory for iid and AR(1) chains", {
  set.seed(10)
  n <- 20000
  iid <- rnorm(n)
  ar <- as.numeric(arima.sim(list(ar = 0.6), n = n))
  e <- effective_size(cbind(iid = iid, ar = ar))
  ## theory: n for iid; n (1 - phi) / (1 + phi) = n / 4 for phi = 0.6.
  ## Relative sampling error of ESS estimates at this n is a few percent.
  expect_lt(abs(e[["iid"]] / n - 1), 0.1)
  expect_lt(abs(e[["ar"]] / (n / 4) - 1), 0.15)
  expect_true(is.na(effective_size(rep(1, 10))))
  if (requireNamespace("coda", quietly = TRUE)) {
    ec <- coda::effectiveSize(ar)
    expect_lt(abs(e[["ar"]] / ec - 1), 0.15)
  }
})

test_that("filter_auc() equals the pairwise definition", {
  set.seed(1)
  s <- round(runif(40), 1)  # ties on purpose
  l <- rbinom(40, 1, 0.5)
  pos <- s[l == 1]
  neg <- s[l == 0]
  ref <- mean(outer(pos, neg, ">")) + 0.5 * mean(outer(pos, neg, "=="))
  expect_equal(filter_auc(l, s), ref)
  expect_true(is.na(filter_auc(c(1, 1), c(0.2, 0.3))))
  expect_error(filter_auc(1, 1:2), "equal length")
})

test_that("overlap_diagnostics() separates overlapping and shifted samples", {
  set.seed(2)
  P <- data.frame(a = rnorm(300), b = rnorm(300))
  same <- data.frame(a = rnorm(300), b = rnorm(300))
  far <- data.frame(a = rnorm(300, 4), b = rnorm(300, -4))
  o1 <- overlap_diagnostics(P, same)
  o2 <- overlap_diagnostics(P, far)
  expect_s3_class(o1, "bsm_overlap")
  expect_gt(o1$OC, 0.8)
  expect_lt(o1$EF, 0.05)
  expect_lt(o1$AUC_PS, 0.6)
  expect_lt(o2$OC, 0.1)
  expect_gt(o2$EF, 0.9)
  expect_gt(o2$AUC_PS, 0.99)
  expect_output(print(o1), "Overlap coefficient")
  expect_error(overlap_diagnostics(P, same[, "a", drop = FALSE]), "missing")
})

test_that("simulate_mixture_data() follows its specification", {
  d <- simulate_mixture_data("binomial", n_prob = 30, n_nonprob = 40, contamination = 0.25,
                             seed = 1)
  expect_equal(nrow(d$prob), 30)
  expect_equal(sum(d$nonprob$contaminated), 10)
  expect_true(all(d$prob$y %in% 0:1))
  expect_equal(d$truth$bad_beta, c(0.2, -1, 0.8))
  d2 <- simulate_mixture_data("binomial", n_prob = 30, n_nonprob = 40, contamination = 0.25,
                              seed = 1)
  expect_identical(d, d2)
  expect_error(simulate_mixture_data(contamination = 2), "contamination")
})
