context("ibdsim")

test_that("ibdsim() returns objects of correct class", {
  x = nuclearPed(1)
  sim = ibdsim(x, N=1, verbose=F, map="uniform.sex.aver", chrom=1:2, model="haldane")
  expect_is(sim, "genomeSimList")
  #expect_is(sim[[1]], "genomeSim")
})
