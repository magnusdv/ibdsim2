context("inbred founders")

quickSim = function(x, sims=1, map="uniform.sex.aver", chromosomes=1, model="haldane")
  ibdsim(x, sims=sims, verbose=F, map=map, chromosomes=chromosomes, model=model)

library(pedtools)

test_that("100% inbred founders are accounted for", {
  x = nuclearPed(1)
  founder_inbreeding(x, 1:2) = 1
  sim = quickSim(x)
  as = alleleSummary(sim[[1]])
  expect_true(all(as[,"3:p"] == 1) && all(as[,"3:m"] == 3))
})

test_that("intermediate inbred founders raise errors", {
  x = nuclearPed(1)
  founder_inbreeding(x, 1) = 0.5
  expect_error(quickSim(x), "Founder inbreeding coefficients other than 0 and 1 are not allowed")
})
