
p1 = function(...) plotSegmentDistribution(..., type = "ibd1", ellipses = F)

pA = function(...) plotSegmentDistribution(..., type = "aut", ellipses = F)

test_that("plotSegmentDistribution(ibd1) catches input errors", {
  x = nuclearPed(2)
  s = ibdsim(x, N=1, map=uniformMap(1), verbose=F, simplify1 = FALSE)
  
  expect_error(p1(s, ids = c(1,5)),
               "Unknown ID label in pedigree 1: 5")
  expect_error(p1(s, ids = list(c(1,5))),
               "Unknown ID label in pedigree 1: 5")
  expect_error(p1(s, ids = list(c(1))),
               "The `ids` entry for pedigree 1 is not a valid pair: 1")
  expect_error(p1(s, ids = list(c(1:3))),
               "The `ids` entry for pedigree 1 is not a valid pair: 1, 2, 3")
})

test_that("plotSegmentDistribution(aut) catches input errors", {
  x = nuclearPed(1)
  s = ibdsim(x, N=1, map=uniformMap(1), verbose=F, simplify1 = FALSE)
  
  expect_error(pA(s, ids = c(1,5)),
               "Unknown ID label in pedigree 1: 5")
  expect_error(pA(s, ids = list(c(1,2))),
               "The `ids` entry for pedigree 1 is not a single individual: 1, 2")
})
