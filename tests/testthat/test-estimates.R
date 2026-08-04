test_that("pairwise estimators check ID length", {
  x = nuclearPed(2)

  expect_error(estimateKinship(x, 3, Nsim = 1), "length 2")
  expect_error(estimateKappa(x, 3, Nsim = 1), "length 2")
  expect_error(estimateIdentity(x, 3, Nsim = 1), "length 2")
  expect_error(
    estimateTwoLocusKappa(x, 3, rho = 0.1, Nsim = 1),
    "length 2")
  expect_error(
    estimateTwoLocusIdentity(x, 3, rho = 0.1, Nsim = 1),
    "length 2")
})


test_that("unlinked estimates use independent simulations", {
  x = nuclearPed(2)

  estimateTwoLocusKappa(x, 3:4, rho = 0.5, Nsim = 5,
                        seed = 1, verbose = FALSE)
  seedAfterTwo = .Random.seed

  estimateKappa(x, 3:4, Nsim = 5, seed = 1, verbose = FALSE)
  seedAfterOne = .Random.seed

  expect_false(identical(seedAfterTwo, seedAfterOne))
})
