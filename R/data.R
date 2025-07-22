#' Legacy version of the `decode19` map
#'
#' A legacy version of the built-in human recombination map, based on
#' Halldorsson et al., 2019. The implementation of this map used in ibdsim2 was
#' updated in v2.3.0, adding the physical endpoint of each chromosome (this was
#' previously lacking), and using a better thinning algorithm to reduce the raw
#' data given by Halldorsson et al. (2019). See also [loadMap()].
#'
#' @format A `genomeMap` object: a list of 23 `chromMap` objects. Each
#'   `chromMap` is a list containing two numeric matrices, named `male` and
#'   `female`, and carries these attributes:
#'
#' * `physStart` – first physical position (Mb) on the chromosome
#' * `physEnd`   – last physical position (Mb) on the chromosome
#' * `physRange` – physical length, `physEnd - physStart` (Mb)
#' * `mapLen`    – length-2 numeric: centiMorgan length of male and female strands
#' * `chrom`     – chromosome label
#' * `Xchrom`    – logical flag used by simulators for X-inheritance
#'
#' @keywords datasets
#'
#' @examples
#' loadMap("decode19")         # new
#' loadMap("legacy_decode19")  # legacy
"legacy_decode19"