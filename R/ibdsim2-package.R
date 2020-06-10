#' ibdsim2: Simulation of chromosomal regions shared by family members
#'
#' Simulation of segments shared identical-by-descent (IBD) by pedigree members.
#' Using sex specific recombination rates along the human genome (Kong et. al
#' (2010) <doi:10.1038/nature09525>), phased chromosomes are simulated for all
#' pedigree members, and patterns of IBD sharing are detected. Additional
#' functions provide further analysis of the simulated genomes.
#'
#' @docType package
#' @import pedtools
#' @importFrom stats rpois runif
#'
#' @name ibdsim2
NULL

#' @useDynLib ibdsim2, .registration = TRUE
#' @importFrom Rcpp sourceCpp evalCpp
NULL