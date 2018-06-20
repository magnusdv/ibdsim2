context("ibdsim")

test_that("ibdsim() returns objects of correct class", {
  x = pedtools::nuclearPed(1)
  sim = ibdsim(x, sims=1, verbose=F, map="uniform.sex.aver", chromosomes=1:2, model="haldane")
  expect_is(sim, "genomeSimList")
  expect_is(sim[[1]], "genomeSim")
  expect_is(sim[[1]][[1]], "chromosomeSim")
  expect_is(sim[[1]][[2]], "chromosomeSim")
})
