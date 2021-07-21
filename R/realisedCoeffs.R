#' Realised relatedness
#'
#' Compute the realised values of various pedigree coefficients, from simulated
#' data. The current implementation covers inbreeding coefficients for single
#' pedigree members, and kinship, kappa and condensed identity coefficients for
#' pairwise relationships.
#'
#' The inbreeding coefficient \eqn{f} of a pedigree member is defined as the
#' probability of autozygosity (homozygous for alleles that are identical by
#' descent) in a random autosomal locus. Equivalently, the inbreeding
#' coefficient is the *expected* autozygous proportion of the autosomal
#' chromosomes.
#'
#' The *realised* inbreeding coefficient \eqn{f_R} in a given individual is the
#' actual fraction of the autosomes covered by autozygous segments. Because of
#' the stochastic nature of meiotic recombination, this may deviate
#' substantially from the pedigree-based expectation.
#'
#' Similarly, the pedigree-based IBD coefficients \eqn{\kappa_0, \kappa_1,
#' \kappa_2} of noninbred pairs of individuals have realised counterparts. For
#' any given pair of individuals we define \eqn{k_i} to be the actual fraction
#' of the autosome where the individuals share exactly \eqn{i} alleles IBD,
#' where \eqn{i = 0,1,2}.
#'
#' Finally, we can do the same thing for each of the nine condensed identity
#' coefficients of Jacquard. For each \eqn{i = 1,...,9} we define \eqn{D_i} the
#' be the fraction of the autosome where a given pair of individuals are in
#' identity state \eqn{i}. This uses the conventional ordering of the nine
#' condensed identity states; see for instance the [`ribd` GitHub
#' page](https://github.com/magnusdv/ribd).

#'
#' @param sims A list of genome simulations, as output by [ibdsim()].
#' @param id,ids A vector with one or two ID labels.
#'
#' @examples
#'
#' # Realised IBD coefficients between full siblings
#' x = nuclearPed(2)
#' s = ibdsim(x, N = 2) # increase N
#' realisedKappa(s, ids = 3:4)
#'
#' ###########
#' 
#' # Realised inbreeding coefficients, child of first cousins
#' x = cousinPed(1, child = TRUE)
#' s = ibdsim(x, N = 2) # increase N
#' realisedInbreeding(s, id = 9)
#'
#' # Same data: realised kinship coefficients between the parents
#' realisedKinship(s, ids = parents(x, 9))
#'
#' ###########
#' 
#' # Realised identity coefficients after full sib mating
#' x = fullSibMating(1)
#' s = ibdsim(x, N = 2) # increase N
#' realisedIdentity(s, ids = 5:6)
#' 
#' @name realised
NULL


#' @rdname realised
#' @importFrom stats sd
#' @export
realisedInbreeding = function(sims, id = NULL) {
  
  # IDs present in sims
  idsims = extractIds(sims)
  
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
      s0 = alleleFlow(s, ids = id, addState = TRUE)
      s = mergeSegments(s0, by = "Aut")
    }

    # Which segments are autozygous
    isAut = as.logical(s[, "Aut"])
    
    # Summary stats on autoz segs
    segLen = s[isAut, 'length']
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
  idsims = extractIds(sims)
  
  if(is.null(ids))
    ids = idsims
  
  if(length(ids) != 2)
    stop2("Argument `ids` must contain exactly two ID labels: ", ids)
  
  if(!all(ids %in% idsims))
    stop2("Target ID not found in segment input:", setdiff(ids, idsims))
  
  if(!is.list(sims))
    sims = list(sims)
  
  resList = vapply(sims, function(s) {
    
    if(!"Sigma" %in% colnames(s)) {
      s0 = alleleFlow(s, ids = ids, addState = TRUE)
      s = mergeSegments(s0, by = "Sigma")
    }
    
    len = s[, 'length']  
    jacq = s[, 'Sigma']  
    
    # Coefficients in the relation phi = weights * jacq
    wei = c(1, 0, .5, 0, .5, 0, .5, .25, 0)
    
    # Weighted sum of segment lengths
    sum(len * wei[jacq])
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
  idsims = extractIds(sims)
  
  if(is.null(ids))
    ids = idsims
  
  if(length(ids) != 2)
    stop2("Argument `ids` must contain exactly two ID labels: ", ids)
  
  if(!all(ids %in% idsims))
    stop2("Target ID not found in segment input:", setdiff(ids, idsims))
  
  if(!is.list(sims))
    sims = list(sims)
  
  resList = lapply(sims, function(s) {
    
    if(!"IBD" %in% colnames(s)) {
      s0 = alleleFlow(s, ids = ids, addState = TRUE)
      s = mergeSegments(s0, by = "IBD")
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
realisedIdentity = function(sims, ids = NULL) {
  
  # IDs present in sims
  idsims = extractIds(sims)
  
  if(is.null(ids))
    ids = idsims
  
  if(length(ids) != 2)
    stop2("Argument `ids` must contain exactly two ID labels: ", ids)
  
  if(!all(ids %in% idsims))
    stop2("Target ID not found in segment input:", setdiff(ids, idsims))
  
  if(!is.list(sims))
    sims = list(sims)
  
  resList = lapply(sims, function(s) {
    
    if(!"Sigma" %in% colnames(s)) {
      s0 = alleleFlow(s, ids = ids, addState = TRUE)
      s = mergeSegments(s0, by = "Sigma")
    }
    
    len = s[, 'length']  
    j = s[, 'Sigma']  
    
    c(D1 = sum(len[j == 1]), 
      D2 = sum(len[j == 2]), 
      D3 = sum(len[j == 3]),
      D4 = sum(len[j == 4]), 
      D5 = sum(len[j == 5]), 
      D6 = sum(len[j == 6]),
      D7 = sum(len[j == 7]), 
      D8 = sum(len[j == 8]), 
      D9 = sum(len[j == 9]),
      nSeg1 = sum(j == 1), 
      nSeg2 = sum(j == 2), 
      nSeg3 = sum(j == 3),
      nSeg4 = sum(j == 4), 
      nSeg5 = sum(j == 5), 
      nSeg6 = sum(j == 6),
      nSeg7 = sum(j == 7), 
      nSeg8 = sum(j == 8), 
      nSeg9 = sum(j == 9))
  })
  
  # Bind row-wise
  resDf = as.data.frame(do.call(rbind, resList))
  
  # Total genome length (assume same for all!)
  L = sum(sims[[1]][, 'length'])
  
  resDf[1:9] = resDf[1:9]/L
  resDf[10:18] = lapply(resDf[10:18], as.integer)
  
  list(perSimulation = resDf, meanCoef = colMeans(resDf[1:9]), stDev = apply(resDf[1:9], 2, sd))
}



#' Extract ID labels from simulation output
#'
#' @param sim Output from [ibdsim()]
#'
#' @return A character vector
#'
#' @examples
#' s = ibdsim(nuclearPed(2), N=1, ids = 3:4)
#' stopifnot(all(extractIds(s) == c("3", "4")))
#' 
#' @export
extractIds = function(sim) {
  clnms = colnames(sim) %||% colnames(sim[[1]])
  if(is.null(clnms))
    stop2("Cannot deduce ID labels from simulation data")
  
  idcols = grep(":p", clnms, fixed = TRUE, value = TRUE)
  substring(idcols, 1, nchar(idcols) - 2)
}
