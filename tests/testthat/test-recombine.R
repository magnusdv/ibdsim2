skip_on_cran()

test_that("recombinations of pure strands", {
  recombine = ibdsim2:::recombine
  s1 = cbind(0,1); s2 = cbind(0,2)
  expect_equal(recombine(s1, s2, numeric()), s1)
  expect_equal(recombine(s1, s2, 0), s2)
  expect_equal(recombine(s1, s2, 1), cbind(c(0,1), c(1,2)))
  expect_equal(recombine(s1, s2, 0:1), cbind(c(0,1), c(2,1)))
  expect_equal(recombine(s1, s2, 1:2), cbind(c(0,1,2), c(1,2,1)))
  expect_equal(recombine(s1, s2, 1:3), cbind(c(0,1,2,3), c(1,2,1,2)))
})

test_that("more recombinations", {
  recombine = ibdsim2:::recombine
  s1 = cbind(c(0,10), c(1,2)); s2 = cbind(0,3)
  expect_equal(recombine(s1, s2, numeric()), s1)
  expect_equal(recombine(s1, s2, 0), s2)
  expect_equal(recombine(s1, s2, 5), cbind(c(0,5), c(1,3)))
  expect_equal(recombine(s1, s2, 10), cbind(c(0,10), c(1,3)))
  expect_equal(recombine(s1, s2, 15), cbind(c(0,10, 15), c(1,2,3)))
})

test_that("back and forth", {
  recombine = ibdsim2:::recombine
  s1 = cbind(c(0,10), c(1,2)); s2 = cbind(c(0,20),c(3,4))
  expect_equal(recombine(s1, s2, c(0.1, 9.9, 10.1)), cbind(c(0, 0.1, 9.9, 10, 10.1, 20),c(1,3,1,2,3,4)))
})
