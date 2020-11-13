
test_that("meiosis results in matrix", {
  pat = cbind(0,1); mat = cbind(0,2)
  map = chromMap(cbind(Mb = c(1,10), cM = c(0,5)))$male
  expect_is(meiosis(list(pat, mat), map), "matrix")
})

test_that("runs of equal alleles are merged", {
  pat = mat = cbind(0, 1)
  map = uniformMap(M = 1)$male
  expect_identical(meiosis(list(pat, mat), map), pat)
})
