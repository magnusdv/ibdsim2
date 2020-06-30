#' ibdsim2: Simulation of chromosomal regions shared by family members
#'
#' Simulation of segments shared identical-by-descent (IBD) by pedigree members.
#' Using sex specific recombination rates along the human genome (Halldorsson et
#' al., 2019), phased chromosomes are simulated for all pedigree members.
#' Additional features include calculation of realised IBD coefficients and IBD
#' segment distribution plots.
#'
#' @docType package
#' @import pedtools
#' @importFrom stats rpois runif
#'
#' @references Halldorsson et al. _Characterizing mutagenic effects of
#'   recombination through a sequence-level genetic map._ Science 363, no. 6425
#'   (2019) \doi{https://doi.org/10.1126/science.aau1043}
#'
#' @name ibdsim2
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp evalCpp
#' @useDynLib ibdsim2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

.onUnload <- function (libpath) {
  library.dynam.unload('ibdsim2', libpath)
}