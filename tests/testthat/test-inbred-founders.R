context("inbred founders")

quickSim = function(x, N=1, map=uniformMap(10), model="haldane")
  ibdsim(x, N=N, verbose=F, map=map, model=model)

test_that("100% inbred founders are accounted for", {
  x = nuclearPed(1)
  founderInbreeding(x, 1:2) = 1
  sim = quickSim(x)[[1]]
  expect_true(all(sim[,"3:p"] == 1) && all(sim[,"3:m"] == 3))
})

test_that("intermediate inbred founders raise errors", {
  x = nuclearPed(1)
  founderInbreeding(x, 1) = 0.5
  expect_error(quickSim(x), "Founder inbreeding coefficients other than 0 and 1 are not allowed")
})
