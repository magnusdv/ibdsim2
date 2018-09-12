context("allele_summary")

library(pedtools)

test_that("allele summary properties of sibs", {
  x = nuclearPed(2)
  sim = ibdsim(x, sims=1, chromosomes=1:2, verbose=F)[[1]]
  summar_all = alleleSummary(sim)
  summar_sibs = alleleSummary(sim, ids=3:4)
  
  expect_is(summar_all, "matrix")
  expect_is(summar_sibs, "matrix")
  
  expect_equal(ncol(summar_all), 4 + 2*pedsize(x))
  expect_equal(ncol(summar_sibs), 4 + 4 + 2)
  
  # parents
  a = alleleSummary(sim, 1:2)
  expect_identical(a[, "IBD"], c(0,0))
  expect_identical(a[, "Sigma"], c(9,9))
})

test_that("allele summary of FSM", {
  x = fullSibMating(1)
  sim = ibdsim(x, sims=1, chromosomes=1, verbose=F, seed=36)[[1]]
  a = alleleSummary(sim, 5:6)
  expect_setequal(a[, 'Sigma'], 1:9)
})
