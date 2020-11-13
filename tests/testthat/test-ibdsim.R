
test_that("ibdsim() returns objects of correct class", {
  x = nuclearPed(1)
  sim = ibdsim(x, N=1, verbose=F, map=uniformMap(1), model="haldane")
  expect_is(sim, "genomeSimList")
})

test_that("Full sib mating yields all 9 jaquard states", {
  x = fullSibMating(1)
  sim = ibdsim(x, N=1, ids = 5:6, map=uniformMap(M=10), verbose=F, seed=1)[[1]]
  expect_setequal(sim[, 'Sigma'], 1:9)
})
