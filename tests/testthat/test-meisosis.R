context("meiosis")

test_that("meiosis results in matrix", {
  pat = cbind(0,1); mat = cbind(0,2)
  map = cbind(Mb=c(1,10), cM=c(0,5))
  expect_is(meiosis(list(pat, mat), map), "matrix")
})
