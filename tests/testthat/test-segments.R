
x = nuclearPed(1)
s = ibdsim(x, N = 1, map = uniformMap(M=1), seed = 123, verbose = F)

test_that("alleleFlow() catches errors", {
  expect_error(alleleFlow(s, ids = 4), "Unknown ID label")
  expect_error(alleleFlow(list(s), ids = 4), 
               "Argument `x` must be a `genomeSim` object. Received: list")
})

test_that("alleleFlow() adds states correctly", {
  ans1 = cbind(chrom=1, startMB=0, endMB=100, startCM=0, endCM=100, `1:p`=1,`1:m`=2)
  ans2 = cbind(chrom=1, startMB=0, endMB=100, startCM=0, endCM=100, `1:p`=1, `1:m`=2, `2:p`=3, `2:m`=4)
  expect_equal(alleleFlow(s, ids = 1, addState = FALSE), ans1) 
  expect_equal(alleleFlow(s, ids = 1, addState = TRUE), cbind(ans1, Aut = 0))
  expect_equal(alleleFlow(s, ids = 1:2, addState = FALSE), ans2) 
  expect_equal(alleleFlow(s, ids = 1:2, addState = TRUE), cbind(ans2, IBD=0, Sigma=9))
  expect_equal(alleleFlow(s, ids = 1:3), s)
})

test_that("zeroIBD() merges adjacent IBD segments before cutoff", {
  sim = cbind(chrom = 1,
              startMB = c(0, 4), endMB = c(4, 8),
              startCM = c(0, 4), endCM = c(4, 8),
              `1:p` = 1, `1:m` = 2,
              `2:p` = 1, `2:m` = 3,
              IBD = 1, Sigma = 8)
  class(sim) = "genomeSim"

  res = zeroIBD(sim, threshold = 8, unit = "mb")
  expect_equal(res, list(zeroprob = 0, stErr = 0))
  expect_error(zeroIBD(sim, threshold = NA), "nonnegative number")
})


test_that("mergeSegments() handles karyogram segment data", {
  # Four-column format produced by prepare_segments()
  segs = data.frame(chrom = 1,
                   startMB = c(0, 1),
                   endMB = c(1, 2),
                   fill = 1)
  res = mergeSegments(segs, by = "fill", checkAdjacency = TRUE)
  expect_equal(nrow(res), 1)
  expect_equal(res$endMB, 2)
})
