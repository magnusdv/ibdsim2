
test_that("ibdsim() returns objects of correct class", {
  x = nuclearPed(1)
  sim1 = ibdsim(x, N=1, verbose=F, map=uniformMap(1), model="haldane", simplify1 = TRUE)
  expect_is(sim1, "genomeSim")
  sim2 = ibdsim(x, N=1, verbose=F, map=uniformMap(1), model="haldane", simplify1 = FALSE)
  expect_is(sim2, "genomeSimList")
})

test_that("Full sib mating yields all 9 jaquard states", {
  x = fullSibMating(1)
  sim = ibdsim(x, ids = 5:6, map=uniformMap(M=10), verbose=F, seed=1)
  expect_setequal(sim[, 'Sigma'], 1:9)
})
