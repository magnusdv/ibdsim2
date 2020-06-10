context("allele_summary")

test_that("allele summary properties of sibs", {
  x = nuclearPed(2)
  sim = ibdsim(x, ids = NULL, N=1, chrom=1:2, verbose=F)[[1]]
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
  sim = ibdsim(x, N=1, ids = 5:6, map=uniformMap(M=10), verbose=F, seed=1)[[1]]
  expect_setequal(sim[, 'Sigma'], 1:9)
})
