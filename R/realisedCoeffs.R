#' Realised relatedness
#'
#' Compute the realised values of various pedigree coefficients, in simulated
#' data. The current implementation covers realised inbreeding coefficients
#' (autozygosity) and realised kappa coefficients.
#'
#' DRAFT: Consider two members A and B of a pedigree P. The *kinship
#' coefficient* between A and B is defined as the probability that a random
#' allele sampled in A is identical by descent (IBD) with an allele sampled in B
#' at the same autosomal locus. If this probability is taken conditional only on
#' the pedigree P, the result is the traditional pedigree-based kinship
#' coefficient.
#'
#' However, because of the discrete nature of meiotic recombination, the actual
#' IBD distribution of individuals with the specified relationship is subject to
#' variation. Hence we may also be interested in the *realised* (or *genomic*)
#' kinship coefficient between A and B, which is the same probability as above,
#' but conditional on the recombination events in the meioses between A and B in
#' P. For example, if the recombination events between A and B happen to result
#' in no segments of IBD sharing, the realised kinship is 0, whatever the
#' pedigree-based kinship may be.
#'
#'
#' @param sims A list of genome simulations, as output by [ibdsim()].
#' @param id A single ID label.
#' @param ids A vector with two ID labels.
#'
#' @examples
#'
#' # Realised IBD coefficients between full siblings
#' x = nuclearPed(2)
#' s = ibdsim(x, N = 10) # increase N
#' realisedKappa(s, ids = 3:4)
#'
#' # Realised inbreeding coefficients, child of first cousins
#' x = cousinPed(1, child = TRUE)
#' s = ibdsim(x, N = 10)
#' realisedAutozygosity(s, id = 9)
#'
#'
#' @name realised
NULL

#' @rdname realised
#' @export
realisedKappa = function(sims, ids = attr(sims, 'ids')) {
  
  if(identical(ids, "leaves"))
    ids = leaves(ped)
  
  if(length(ids) != 2)
    stop2("Argument `ids` must contain exactly two ID labels: ", ids)

  L = attr(sims, 'genome_length_Mb')
  
  segs = vapply(sims, function(s) {
    a = segmentSummary(s, ids, addState = TRUE)
    
    # Merge adjacent segments with equal IBD state
    a2 = mergeAdjacent(a, vec = a[, 'IBD'])

    len = a2[, 'length']  
    ibd = a2[, 'IBD']  
      
    c(ibd0 = sum(len[ibd == 0]), 
      ibd1 = sum(len[ibd == 1]), 
      ibd2 = sum(len[ibd == 2]), 
      Nseg1 = sum(ibd == 1), 
      Nseg2 = sum(ibd == 2), 
      Nseg = sum(ibd > 0))
  }, FUN.VALUE = numeric(6))
  
  kappaRealised = segs[1:3, , drop = FALSE]/L
  list(kappaRealised = kappaRealised, 
       nSegments = segs[4:6, , drop = FALSE], 
       kappaHat = rowMeans(kappaRealised),
       genomeLength = L)
}


#' @rdname realised
#' @export
realisedAutozygosity = function(sims, id = attr(sims, 'ids')) {
  # Pedigree
  ped = attr(sims, 'pedigree')
  
  if(identical(id, "leaf"))
    id = leaves(ped)
  
  if(length(id) != 1)
    stop2("Argument `id` must contain a single ID label: ", id)
  
  # Theoretical inbreeding coefficient
  inb = ribd::inbreeding(ped)[[internalID(ped, id)]] # `[[` to remove name
  
  # Total genome length
  genomeL = attr(sims, "genome_length_Mb")
  
  # Summarise each simulation (previously used vapply(); now lapply + do.call(rbind))
  summList = lapply(sims, function(s) {
    a = segmentSummary(s, ids = id, addState = TRUE)
    
    # Merge adjacent segments with equal IBD state
    a2 = mergeAdjacent(a, vec = as.logical(a[, 'Aut']))
    
    aut = as.logical(a2[, "Aut"])
    segLengths = a2[aut, 'length']
    segCount = length(segLengths)
    totLength = sum(segLengths)
    meanLen = ifelse(segCount == 0, 0, totLength/segCount)
    
    c(segCount = segCount, meanLength = meanLen, 
      totLength = totLength, 
      longest = max(segLengths), shortest = min(segLengths),
      fraction = totLength/genomeL)
  })
  
  summDF = as.data.frame(do.call(rbind, summList), make.names = FALSE)
  
  cbind(summDF, expected = inb, genomeLength = genomeL)
}
