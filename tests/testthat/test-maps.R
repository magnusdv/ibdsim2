context("plot segment distributions")

test_that("uniformMap() returns a chromMap of correct length", {
  m = uniformMap(Mb = 1, cM = 2:3)
  expect_s3_class(m, "chromMap")
  expect_true(isChromMap(m))
  expect_true(is.data.frame(m$male))
  expect_true(is.data.frame(m$female))
  expect_equal(chromLen(m), 1)
  expect_equal(chromLen(m, "cM"), c(male=2, female=3))
  expect_equal(chromLen(m, "cM", "male"), 2)
  expect_equal(chromLen(m, "cM", "female"), 3)
})

test_that("genomeMap() returns a genomeMap of correct length", {
  m = uniformMap(Mb = 1, cM = 2:3)
  g = genomeMap(m)
  expect_s3_class(g, "genomeMap")
  expect_true(isGenomeMap(g))
  expect_equal(g[[1]], m)
  expect_equal(genomeLen(g), 1)
  expect_equal(genomeLen(g, "cM"), c(male=2, female=3))
  expect_equal(genomeLen(g, "cM", "male"), 2)
  expect_equal(genomeLen(g, "cM", "female"), 3)
})

test_that("loadMap() catches errors", {
  expect_error(loadMap(""), "Unknown map")
  expect_error(loadMap(1), "Argument `map` must be a character of length 1")
  expect_error(loadMap(letters[1:2]), "Argument `map` must be a character of length 1")
  expect_error(loadMap(detailed = 1), "Argument `detailed` must be either TRUE or FALSE")
  expect_error(loadMap(detailed = NULL), "Argument `detailed` must be either TRUE or FALSE")
  expect_error(loadMap(sex = NA), "Argument `sexSpecific` must be either TRUE or FALSE")
  expect_error(loadMap(sex = c(T,T)), "Argument `sexSpecific` must be either TRUE or FALSE")
  expect_error(loadMap(chrom = 30), "Index out of range")
})


