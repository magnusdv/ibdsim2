xprof = nuclearPed(1) |>
  addMarker(alleles = 1:100, chrom = 1, posMb = 50) |>
  addMarker(alleles = 1:100, chrom = 1, posMb = 25)

sprof = cbind(
  chrom = 1,
  startMB = c(0, 50),
  endMB = c(50, 100),
  startCM = c(0, 50),
  endCM = c(50, 100),
  `1:p` = 1,
  `1:m` = 1,
  `2:p` = c(1, 2),
  `2:m` = c(1, 2))
class(sprof) = "genomeSim"


test_that("profileSimIBD() follows the IBD pattern", {
  set.seed(123)
  als = sample.int(100, size = 2, replace = TRUE,
                   prob = rep(0.01, 100))
  a1 = als[1]
  a2 = als[2]
  y = profileSimIBD(xprof, sprof, ids = 1:2, seed = 123, verbose = FALSE)

  # Marker exactly at 50 Mb uses the second segment
  expect_equal(y$MARKERS[[1]][1, ], c(a1, a1))
  expect_equal(y$MARKERS[[1]][2, ], c(a2, a2))

  # Both founders share code 1 in the first segment
  expect_equal(length(unique(as.vector(y$MARKERS[[2]][1:2, ]))), 1)

  # Individuals outside ids are cleared
  expect_equal(y$MARKERS[[1]][3, ], c(0, 0))
})


test_that("profileSimIBD() handles lists and marker selection", {
  y = profileSimIBD(xprof, list(sprof, sprof), ids = 1:2,
                    markers = 1, seed = 1, verbose = FALSE)

  expect_length(y, 2)
  expect_equal(unname(sapply(y, nMarkers)), c(1, 1))
})


test_that("profileSimIBD() reports missing chromosomes", {
  bad = sprof
  bad[, "chrom"] = 2

  expect_error(
    profileSimIBD(xprof, bad, ids = 1:2, verbose = FALSE),
    "Chromosome missing from `ibdpattern`: 1")
})


test_that("profileSimIBD() handles X males", {
  x = nuclearPed(1, sex = 1) |>
    addMarker(alleles = 1:4, chrom = "X", posMb = 10)

  sim = cbind(chrom = 23, startMB = 0, endMB = 100,
              startCM = 0, endCM = 100,
              `3:p` = 0, `3:m` = 1)
  class(sim) = "genomeSim"

  y = profileSimIBD(x, sim, seed = 1, verbose = FALSE)
  expect_equal(y$MARKERS[[1]][3, 1], y$MARKERS[[1]][3, 2])
})

