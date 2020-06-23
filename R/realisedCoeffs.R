#' Realised relatedness
#'
#' Compute the realised values of various pedigree coefficients, from simulated
#' data. The current implementation covers realised inbreeding coefficients for
#' single pedigree members, and realised kinship and kappa coefficients for
#' pairwise relationships.
#'
#' The inbreeding coefficient \eqn{f} of a pedigree member is defined as the
#' probability of autozygosity (homozygous for alleles that are identical by
#' descent) in a random autosomal locus. Equivalently, the inbreeding
#' coefficient is the *expected* autozygous proportion of the autosomal
#' chromosomes.
#'
#' The *realised* inbreeding coefficient \eqn{f_real} in a given individual is
#' the actual fraction of the autosomes covered by autozygous segments. Because
#' of the stochastic nature of meiotic recombination, this may deviate
#' substantially from the pedigree-based expectation.
#'
#' Similarly, the pedigree-based IBD coefficients \eqn{\kappa = (\kappa_0,
#' \kappa_1, \kappa_2)} have realised counterparts, when looking at a specific
#' pair of individuals:
#'
#' * \eqn{k_0}: The actual fraction of the autosome where the individuals share
#' 0 alleles IBD
#'
#' * \eqn{k_1}: The actual fraction of the autosome where the individuals share
#' 1 alleles IBD
#'
#' * \eqn{k_2}: The actual fraction of the autosome where the individuals share
#' 2 alleles IBD
#'
#' @param sims A list of genome simulations, as output by [ibdsim()].
#' @param id,ids A vector with one or two ID labels.
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
#' s = ibdsim(x, N = 10) # increase N
#' realisedInbreeding(s, id = 9)
#'
#' # Same data: realised kinship coefficients between the parents
#'  realisedKinship(s, ids = parents(x, 9))
#'
#' @name realised
NULL


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


#' @rdname realised
#' @importFrom stats sd
#' @export
realisedKinship = function(sims, ids = NULL) {
  
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
  
  resList = vapply(sims, function(s) {
    
    if(length(idsims) > 2 || !"IBD" %in% colnames(s)) {
      s0 = segmentSummary(s, ids = ids, addState = TRUE)
      s = mergeAdjacent(s0, vec = "IBD")
    }
    
    len = s[, 'length']  
    ibd = s[, 'IBD']  
    
    ibd1 = sum(len[ibd == 1])
    ibd2 = sum(len[ibd == 2])
    
    ibd1/4 + ibd2/2
  }, FUN.VALUE = numeric(1))
  
  # Total genome length (assume same for all!)
  L = sum(sims[[1]][, 'length'])
  
  resVec = unlist(resList)/L
  
  list(perSimulation = resVec, meanCoef = mean(resVec), stDev = sd(resVec))
}


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


# Utility for deducing ID labels present in an output of `ibdsim()` simulation
extractIdsFromSegmentSummary = function(x) {
  clnms = colnames(x) %||% colnames(x[[1]])
  if(is.null(clnms))
    stop2("Cannot deduce ID labels from segment summary")
  
  idcols = grep(":p", clnms, fixed = TRUE, value = TRUE)
  substring(idcols, 1, nchar(idcols) - 2)
}
