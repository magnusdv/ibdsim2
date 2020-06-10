context("plot IBD1")

pp = function(...) plotIBD1(..., ellipses = F)

test_that("plot_ibd1() catches input errors", {
  x = nuclearPed(2)
  s = ibdsim(x, sims=10, chromosomes=21, verbose=F)
  expect_error(pp(s, pairs = list(PO = c(1,5))),
               "Unknown ID label in pedigree PO: 5")
  expect_error(pp(s, pairs = list(PO = c(1))),
               "The `pairs` entry for pedigree PO is not a valid pair: 1")
  expect_error(pp(s, pairs = list(PO = c(1:3))),
               "The `pairs` entry for pedigree PO is not a valid pair: 1, 2, 3")
  expect_error(pp(s, pairs = c(1,5)),
               "Unknown ID label in pedigree 1: 5")
  expect_error(pp(s, pairs = list(c(1,5))),
               "Unknown ID label in pedigree 1: 5")
  expect_error(pp(s, pairs = list(c(1))),
               "The `pairs` entry for pedigree 1 is not a valid pair: 1")
  expect_error(pp(s, pairs = list(c(1:3))),
               "The `pairs` entry for pedigree 1 is not a valid pair: 1, 2, 3")
})
