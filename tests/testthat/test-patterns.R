
# Setup example
x = fullSibMating(1)
s = ibdsim(x, N = 1, map = uniformMap(M=1), seed = 1234, verbose = F)[[1]]


test_that("findPattern() catches errors", {
  expect_error(findPattern(s, list(FOO = 1:3)), "should be one of")
  expect_error(findPattern(s, list(aut = 3:4, non = 3)), " both carrier and noncarrier: 3")
  expect_error(findPattern(s, list(aut = 3:4, het = 3)), " both autozygous and heterozygous: 3")
  expect_error(findPattern(s, list(aut = "foo")), "Unknown ID label: foo")
  
})

test_that("findPattern() finds segments correctly", {
  expect_equal(nrow(findPattern(s, list(aut = 6))), 2)
  expect_equal(nrow(findPattern(s, list(aut = 5:6))), 1)
  expect_equal(nrow(findPattern(s, list(aut = 1:2))), 0)
  expect_equal(nrow(findPattern(s, list(aut = 6, non = 1))), 2)
  expect_equal(nrow(findPattern(s, list(aut = 6, non = 2))), 0)
  
  expect_equal(nrow(findPattern(s, list(aut = 6, het = 5))), 2)
  expect_equal(nrow(findPattern(s, list(carrier = c(1,6), non = 5))), 2)
})

test_that("findPattern() merges segments correctly", {
  expect_equal(nrow(findPattern(s, list(het = 3), merge = FALSE)), 3)
  expect_equal(nrow(findPattern(s, list(het = 3), merge = TRUE)), 1)
  
})
