#' Allele sharing summary
#'
#' This function facilitates downstream analysis of simulations produced by \code{\link{ibdsim}}.
#' It summarises a single genome simulation by describing the allele flow through the pedigree.
#' 
#' @param x An object of class \code{genomeSim}, i.e. a list of simulated chromosomes. 
#' Each chromosome is a list, with one entry for each individual. Each of these 
#' entries is a list of two matrices (one for each strand). The matrices have 2 
#' columns (start position; allele) and one row for each segment unbroken by recombination.
#' @param ids A vector of numerical IDs. If missing, all individuals are included.
#' @param ibd.status TRUE or FALSES. This parameter is meaningful only if 
#' \code{length(ids)==2}. If TRUE the IBD status (number of alleles shared IBD, either 
#' 0, 1 or 2) of each segment is computed, as well as the breakdown of their parental origin.

#'
#' @return A numerical matrix. Each row corresponds to a chromosomal segment. 
#' The first 4 columns describe the segment (chromosome, start, end, length), 
#' and are followed by two columns (paternal allele, maternal allele) for each 
#' of the selected individuals. If \code{ibd.status=TRUE} five more columns are
#' added: \code{ibd}, \code{ibd_pp}, \code{ibd_pm}, \code{ibd_mp} and \code{ibd_mm}. 
#' The first of these indicate the IBD status (0, 1 or 2) in the segment, 
#' while the latter 4 give the parental breakdown of this number. For instance,
#' \code{ibd_pm} is 1 if the _p_aternal allele of the first individual is IBD 
#' with the _m_aternal allele of the second individual, and 0 otherwise.
#'
#' @examples
#' ### Sibling simulation (3 sims of chromosomes 1 and 2)
#' x = pedtools::nuclearPed(2)
#' sim = ibdsim(x, sims=3, chromosomes=1:2)
#' 
#' alleleSummary(sim[[1]]) # First sim, summary of all individuals
#' alleleSummary(sim[[1]], ids=3:4) # Summary of the siblings
#' alleleSummary(sim[[1]], ids=3:4, ibd.status=TRUE) # IBD breakdown of the siblings
#' 
#' # Trivial example: Summary of the father.
#' # Being the first founder, his alleles are denoted 1 and 2.
#' alleleSummary(sim[[1]], ids=1) 
#' 
#' @importFrom pedtools internalID
#' @export
alleleSummary = function(x, ids, ibd.status=FALSE) {
  assert_that(inherits(x, "genomeSim"))
  
  ped = attr(x, "pedigree")
  if (missing(ids)) 
    ids = ped$LABELS
  ids = internalID(ped, ids) # ad hoc. TODO
  
  if(ibd.status && length(ids)!=2)
    stop("Parameter 'ibd.status' is meaningful only if length(ids)==2.")
  
  allele.colnames = paste0(rep(ids, each = 2), c("p", "m"))
print(ids)
  each.chrom = lapply(x, function(y) {
    haplos = unlist(y[ids], recursive = FALSE)
    breaks = unlist(lapply(haplos, function(m) m[-1, 1]))
    breaks = c(0, .sortDouble(breaks[!duplicated(breaks)]))
    alleles = vapply(haplos, pos2allele, posvec = breaks, FUN.VALUE = breaks)
    if (length(breaks) == 1) # Ad hoc fix, since vapply simplifies if FUN.VALUE has length 1
      dim(alleles) = c(1, length(allele.colnames))
    colnames(alleles) = allele.colnames

    stops = c(breaks[-1], attr(y, "length_Mb"))
    chrom = rep.int(attr(y, "chromosome"), length(breaks))
    cbind(chrom = chrom, start = breaks, end = stops, length = stops - breaks, alleles)
  })
  res = do.call(rbind, each.chrom)
  if (ibd.status) {
    ibd_pp = res[, 5] == res[, 7]
    ibd_pm = res[, 5] == res[, 8]
    ibd_mp = res[, 6] == res[, 7]
    ibd_mm = res[, 6] == res[, 8]
    ibd = (ibd_pp | ibd_pm) + (ibd_mp | ibd_mm)

    res = cbind(res, ibd = ibd, ibd_pp = ibd_pp, ibd_pm = ibd_pm, ibd_mp = ibd_mp, ibd_mm = ibd_mm)
  }
  res
}
