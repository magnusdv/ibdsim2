# library call causes annoying startup messages
# library(pedtools, quietly = TRUE)

nuclearPed = pedtools::nuclearPed
addSon = pedtools::addSon
addDaughter = pedtools::addDaughter
cousinPed = pedtools::cousinPed
halfSibPed = pedtools::halfSibPed
linearPed = pedtools::linearPed
halfCousinPed = pedtools::halfCousinPed

BUILTIN_PEDS = list(        
  "Trio" = nuclearPed(1),
  "Siblings" = nuclearPed(2),
  "Sibship of 3" = nuclearPed(3, sex = c(1,2,1)),
  "Half-sibs, maternal" = halfSibPed(1, 1, type = "maternal"),
  "Half-sibs, paternal" = halfSibPed(1, 1),
  "Avuncular, maternal" = nuclearPed(2, sex = 1:2) |> addSon(4, verbose = FALSE),
  "Avuncular, paternal" = nuclearPed(2) |> addSon(4, verbose = FALSE),
  "Grandparent (female line)" = linearPed(2, sex = 2),
  "Grandparent (male line)" = linearPed(2),
  "Great grandp (female line)" = linearPed(3, sex = 2),
  "Great grandp (male line)" = linearPed(3),
  "1st cousins" = cousinPed(1, symmetric = TRUE),
  "1st cousins + child" = cousinPed(1, symmetric = TRUE, child = TRUE),
  "2nd cousins" = cousinPed(2, symmetric = TRUE),
  "2nd cousins + child" = cousinPed(2, symmetric = TRUE, child = TRUE),
  "Half 1st cousins" = halfCousinPed(1, symmetric = TRUE),
  "Half 1st cousins + child" = halfCousinPed(1, symmetric = TRUE, child = TRUE),
  "Half 2nd cousins" = halfCousinPed(2, symmetric = TRUE),
  "Half 2nd cousins + child" = halfCousinPed(2, symmetric = TRUE, child = TRUE),
  "3/4-siblings" = nuclearPed(2) |> addSon(c(3,5), verbose = FALSE) |> addSon(4:5),
  "3/4-siblings + child" = nuclearPed(2) |> addSon(c(3,5), verbose = FALSE) |> addDaughter(4:5) |> addSon(6:7),
  "Dbl 1st cousins" = pedtools::doubleFirstCousins(),
  "Dbl 1st cousins + child" = pedtools::doubleCousins(1, 1, child = TRUE),
  "Quad half 1st cousins" = pedtools::quadHalfFirstCousins(),
  "Father-daughter incest" = nuclearPed(1, sex = 2) |> addSon(c(1,3)),
  "Mother-son incest" = nuclearPed(1) |> addSon(2:3),
  "Full-sib incest" = nuclearPed(2, sex = 1:2) |> addSon(3:4),
  "Half-sib incest" = halfSibPed(sex2 = 2) |> addSon(4:5),
  "Grandfather incest" = linearPed(2, sex = 1:2) |> addSon(c(1,5)),
  "Grandmother incest" = linearPed(2, sex = 2:1) |> addSon(c(2,5))
)

# Default individuals to be indicated when a pedigree is loaded
DEFAULT_IDS = lapply(BUILTIN_PEDS, pedtools::leaves)

# A few special cases
DEFAULT_IDS[["Grandparent (female line)"]] = c(2,5)
DEFAULT_IDS[["Grandparent (male line)"]] = c(1,5)
DEFAULT_IDS[["Great grandp (female line)"]] = c(2,7)
DEFAULT_IDS[["Great grandp (male line)"]] = c(1,7)
