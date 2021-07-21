
x = nuclearPed(1)
s = ibdsim(x, N = 1, map = uniformMap(M=1), seed = 123, verbose = F)[[1]]

test_that("alleleFlow() catches errors", {
  expect_error(alleleFlow(s, ids = 4), 
               "Unknown ID label")
  expect_error(alleleFlow(list(s), ids = 4), 
               "Argument `x` must be a `genomeSim` object. Received: list")
})

test_that("alleleFlow() adds states correctly", {
  ans1 = cbind(chrom=1, start=0, end=100, length=100, `1:p`=1,`1:m`=2)
  ans2 = cbind(chrom=1, start=0, end=100, length=100, `1:p`=1, `1:m`=2, `2:p`=3, `2:m`=4)
  expect_equal(alleleFlow(s, ids = 1, addState = FALSE), ans1) 
  expect_equal(alleleFlow(s, ids = 1, addState = TRUE), cbind(ans1, Aut = 0))
  expect_equal(alleleFlow(s, ids = 1:2, addState = FALSE), ans2) 
  expect_equal(alleleFlow(s, ids = 1:2, addState = TRUE), cbind(ans2, IBD=0, Sigma=9))
  expect_equal(alleleFlow(s, ids = 1:3), s)
})
