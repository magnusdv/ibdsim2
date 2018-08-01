context("allele_summary")

test_that("allele summary properties of sibs", {
  x = pedtools::nuclearPed(2)
  sim1 = ibdsim(x, sims=1, chromosomes=1:2, verbose=F)[[1]]
  summar_all = alleleSummary(sim1)
  summar_sibs = alleleSummary(sim1, ids=3:4, ibd=T)
  
  expect_is(summar_all, "matrix")
  expect_is(summar_sibs, "matrix")
  
  expect_equal(ncol(summar_all), 4 + 2*pedsize(x))
  expect_equal(ncol(summar_sibs), 4 + 4 + 5)
  
  expect_true(all(summar_sibs[,'ibd'] == rowSums(summar_sibs[,10:13])))
})
