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
#' realisedInbreeding(s, id = 9)
#'
#'
#' @name realised
NULL

#' @rdname realised
#' @importFrom stats sd
#' @export
realisedKappa = function(sims, ids = NULL) {
  
  # IDs present in sims
  idsims = extractIdsFromSegmentSummary(sims)
  
  if(is.null(ids))
    ids = idsims
  
  if(length(ids) != 2)
    stop2("Argument `ids` must contain exactly two ID labels: ", ids)
  
  if(!all(ids %in% idsims))
    stop2("Target ID not found in segment input:", setdiff(ids, idsims))
  
  if(!is.list(sims))
    sims = list(sims)
  
  resList = lapply(sims, function(s) {
    
    if(length(idsims) > 2 || !"IBD" %in% colnames(s)) {
      s0 = segmentSummary(s, ids = ids, addState = TRUE)
      s = mergeAdjacent(s0, vec = "IBD")
    }
    
    len = s[, 'length']  
    ibd = s[, 'IBD']  
      
    c(ibd0 = sum(len[ibd == 0]), ibd1 = sum(len[ibd == 1]), ibd2 = sum(len[ibd == 2]),
      nSeg0 = sum(ibd == 0), nSeg1 = sum(ibd == 1), nSeg2 = sum(ibd == 2))
  })
  
  # Bind row-wise
  resMat = do.call(rbind, resList)
  
  # Total genome length (assume same for all!)
  L = sum(sims[[1]][, 'length'])
  
  resDf = data.frame(k0 = resMat[,1]/L, 
                     k1 = resMat[,2]/L, 
                     k2 = resMat[,3]/L, 
                     nSeg0 = as.integer(resMat[,4]), 
                     nSeg1 = as.integer(resMat[,5]), 
                     nSeg2 = as.integer(resMat[,6]))
  
  list(perSimulation = resDf, meanCoef = colMeans(resDf[, 1:3]), stDev = apply(resDf[,1:3], 2, sd))
}


#' @rdname realised
#' @importFrom stats sd
#' @export
realisedInbreeding = function(sims, id = NULL) {
  
  # IDs present in sims
  idsims = extractIdsFromSegmentSummary(sims)
  
  # Target ID
  if(is.null(id))
    id = idsims
  if(length(id) != 1)
    stop2("Argument `id` must contain a single ID label: ", id)
  if(!id %in% idsims)
    stop2("Target ID not found in the IBD segment input")
  
  # Ensure sims is a list
  if(!is.null(dim(sims)))
    sims = list(sims)
  
  # Summarise each simulation
  resList = lapply(sims, function(s) {
    
    if(length(idsims) > 1 || !"Aut" %in% colnames(s)) {
      s0 = segmentSummary(s, ids = id, addState = TRUE)
      s = mergeAdjacent(s0, vec = "Aut")
    }

    # Which segments are autozygous
    aut = as.logical(s[, "Aut"])
    
    # Summary stats on autoz segs
    segLen = s[aut, 'length']
    nSeg = length(segLen)
    totLen = sum(segLen)
    
    c(nSeg = nSeg, 
      meanLen = ifelse(nSeg == 0, 0, totLen/nSeg), 
      totLen = totLen, 
      maxLen = ifelse(nSeg == 0, 0, max(segLen)), 
      minLen = ifelse(nSeg == 0, 0, min(segLen)))
  })
  
  # Bind row-wise
  resDf = as.data.frame(do.call(rbind, resList), make.names = FALSE)
  
  # Total genome length (assume same for all!)
  L = sum(sims[[1]][, 'length'])
  
  # Add realised f
  resDf$fReal = fr = resDf$totLen/L
  
  list(perSimulation = resDf, meanCoef = mean(fr), stDev = sd(fr))
}

# Utility for deducing ID labels present in an output of `ibdsim()` simulation
extractIdsFromSegmentSummary = function(x) {
  clnms = colnames(x) %||% colnames(x[[1]])
  if(is.null(clnms))
    stop2("Cannot deduce ID labels from segment summary")
  
  idcols = grep(":p", clnms, fixed = TRUE, value = TRUE)
  substring(idcols, 1, nchar(idcols) - 2)
}
