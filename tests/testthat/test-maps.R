
test_that("uniformMap() returns a chromMap of correct length", {
  m = uniformMap(Mb = 1, cM = 2:3)
  expect_s3_class(m, "chromMap")
  expect_true(isChromMap(m))
  expect_true(is.data.frame(m$male))
  expect_true(is.data.frame(m$female))
  expect_equal(physRange(m), 1)
  expect_equal(mapLen(m,), c(male=2, female=3))
  expect_equal(mapLen(m, "male"), 2)
  expect_equal(mapLen(m, "female"), 3)
})

test_that("genomeMap() returns a genomeMap of correct length", {
  m = uniformMap(Mb = 1, cM = 2:3)
  g = genomeMap(m)
  expect_s3_class(g, "genomeMap")
  expect_true(isGenomeMap(g))
  expect_equal(g[[1]], m)
  expect_equal(physRange(g), 1)
  expect_equal(mapLen(g), c(male=2, female=3))
  expect_equal(mapLen(g, "male"), 2)
  expect_equal(mapLen(g, "female"), 3)
})

test_that("loadMap() catches errors", {
  expect_error(loadMap(""), "Unknown map")
  expect_error(loadMap(1), "Argument `map` must be a character of length 1")
  expect_error(loadMap(letters[1:2]), "Argument `map` must be a character of length 1")
  expect_error(loadMap(uniform = 1), "Argument `uniform` must be either TRUE or FALSE")
  expect_error(loadMap(uniform = NULL), "Argument `uniform` must be either TRUE or FALSE")
  expect_error(loadMap(sex = NA), "Argument `sexAverage` must be either TRUE or FALSE")
  expect_error(loadMap(sex = c(T,T)), "Argument `sexAverage` must be either TRUE or FALSE")
  expect_error(loadMap(chrom = 30), "Unknown chromosome name")
  expect_error(loadMap(chrom = NA), "Unknown chromosome name")
  expect_error(loadMap(chrom = c(1,1)), "Duplicated chromosome")
  expect_error(loadMap(chrom = c(23,"X")), "Duplicated chromosome")
})

test_that("loadMap() options work", {
  m = loadMap(chrom = 13)
  mu = loadMap(chrom = 13, uniform = T)
  ms = loadMap(chrom = 13, sexAver = T)
  mus = loadMap(chrom = 13, uniform = T, sexAver = T)
  
  expect_length(m, 1)
  expect_equal(physRange(mu), physRange(m))
  expect_equal(physRange(ms), physRange(m))
  expect_equal(physRange(mus), physRange(m))
  expect_equal(mapLen(mu), mapLen(m))
  expect_equal(mapLen(ms), mapLen(mus))
  expect_equal(mapLen(mus, "male"), mean(mapLen(m)))
  expect_equal(nrow(mu[[1]]$male), 2)
  expect_equal(nrow(mus[[1]]$male), 2)
  expect_equal(nrow(m[[1]]$male), nrow(ms[[1]]$male))
})

test_that("customMap() results in a genomeMap", {
  expect_is(customMap(data.frame(CHROM = 30, MB = 0, CM = 0, blabla = NA)), "genomeMap")
})

test_that("customMap() catches errors", {
  expect_error(customMap(list(chrom = 30, foo = 0, bar = 0)), 
               "Argument `x` must be a data frame or matrix.")
  expect_error(customMap(data.frame(foo_chrom = 30, foo = 0, bar = 0)), 
               '`x` must have a column named "chrom"')
  expect_error(customMap(data.frame(chrom = 30, foo = 0, bar = 0)), 
               '`x` must have a column named "mb"')
  expect_error(customMap(data.frame(chrom = 30, mb = 0, male = 0)), 
               '`x` must either have a column named "cm", or two columns named "male" and "female".')
})

test_that("customMap() assigns male/female columns correctly", {
  m1 = customMap(data.frame(chrom = 1, mb = 0:1, male = c(0,2), female = c(0,3)))
  expect_equal(mapLen(m1), c(male = 2, female = 3))
  
  m2 = customMap(data.frame(chrom = 1, mb = 0:1, female = c(0,3), male = c(0,2)))
  expect_equal(mapLen(m2), c(male = 2, female = 3))
})

test_that("uniformMap() catches ambiguous and invalid lengths", {
  expect_error(uniformMap(cM = 2:3), "`Mb` must be supplied")
  expect_error(uniformMap(Mb = 0, cM = 1), "Genetic map length")

  m = uniformMap(Mb = 1, cM = c(2, 2))
  expect_equal(mapLen(m), c(male = 2, female = 2))
})


test_that("convertPos() handles vectors and X maps", {
  m = uniformMap(Mb = 100, cM = 100)
  g = genomeMap(m)

  expect_equal(convertPos(chrom = 1, Mb = c(0, 50, 100), map = g),
               c(0, 50, 100))

  mx = uniformMap(Mb = 100, cM = 100, chrom = "X")
  gx = genomeMap(mx)

  expect_equal(convertPos(Mb = c(0, 50, 100), map = mx),
               c(0, 50, 100))
  expect_equal(convertPos(chrom = 23, Mb = 50, map = gx), 50)
  expect_error(convertPos(Mb = 50, map = mx, sex = "male"),
               "no male component")
})
