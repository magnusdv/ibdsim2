#' Find specific IBD patterns
#'
#' Find segments satisfying a particular pattern of IBD sharing, in a list of
#' IBD simulations.
#'
#' For each simulation, this function extracts the subset of rows satisfying the
#' allele sharing specified by `pattern`. That is, segments where, for some allele,
#'
#' * all of `pattern$autozygous` are autozygous
#' 
#' * all of `pattern$heterozygous` have exactly one copy
#' 
#' * all of `pattern$carriers` have at least one copy
#' 
#' * none of `pattern$noncarriers` carry the allele.
#' 
#'
#' @param sims A `genomeSim` object, or a list of such. Typically made by
#'   [ibdsim()].
#' @param pattern A named list of vectors containing ID labels. Allowed names
#'   are `autozygous`, `heterozygous`, `carriers`, `noncarriers`.
#' @param merge A logical, indicating if adjacent segments should be merged.
#'   Default: TRUE.
#' @param cutoff A non-negative number. Segments shorter than this are excluded
#'   from the output. Default: 0.
#'
#' @return A matrix (if `sims` is a single `genomeSim` object), or a list of
#'   matrices.
#'
#' @seealso [segmentStats()]
#' 
#' @examples
#' x = nuclearPed(3)
#' s = ibdsim(x, N = 1, map = uniformMap(M = 1), seed = 1729)
#' s1 = s[[1]]
#'
#' # Segments where some allele is shared by 3 and 4, but not 5
#' pattern = list(carriers = 3:4, noncarriers = 5)
#' findPattern(s1, pattern)
#'
#' # Exclude segments less than 7 cM
#' findPattern(s1, pattern, cutoff = 7)
#'
#' # Visual confirmation:
#' haploDraw(x, s1, margin = c(5,3,3,3))
#'
#' @export
findPattern = function(sims, pattern, merge = TRUE, cutoff = 0) {
  
  if(single <- is.matrix(sims))
    sims = list(sims)
  
  names(pattern) = match.arg(names(pattern), 
                             c("autozygous", "heterozygous", "carriers", "noncarriers"), 
                             several.ok = TRUE)
  
  pattern = lapply(pattern, as.character)
  
  aut = pattern$autozygous
  het = pattern$heterozygous
  car = setdiff(pattern$carriers, c(aut, het))
  non = pattern$noncarriers
  
  # All individuals involved
  allids = c(aut, het, car, non)
  
  # All carriers
  allcarr = c(aut, het, car)
  
  # If no carriers, return. TODO: Why return all of sims???
  if(!length(allcarr)) 
    return(sims)

  if(length(err <- intersect(non, allcarr)))
    stop2("Individuals indicated as both carrier and noncarrier: ", err)
  
  if(length(err <- intersect(aut, het)))
    stop2("Individuals indicated as both autozygous and heterozygous: ", err)
  
  # First carrier (which group is irrelevant)
  id1 = allcarr[1]
  
  # List of allele columns of each individual
  cols = lapply(allids, paste0, c(":p", ":m"))
  names(cols) = allids
  
  # Utils extracting alleles
  pat = function(sim, id) sim[, cols[[id]][1]]
  mat = function(sim, id) sim[, cols[[id]][2]]
  
  # Loop over simulations
  res = lapply(sims, function(sim) {
    s = alleleFlow(sim, allids, addState = FALSE)
    
    # The possible alleles to be shared
    a1 = pat(s, id1)
    a2 = mat(s, id1)
    
    # Initialise logical vectors indicating matching rows
    eq1 = eq2 = rep(TRUE, length = nrow(s))
    
    ### Find rows satisfying pattern
    
    # Autozygous carriers
    for(id in aut) {
      p = pat(s, id)
      isAut = p == mat(s, id)
      eq1 = eq1 & isAut & a1 == p
      eq2 = eq2 & isAut & a2 == p
    }
    
    # Heterozygous carriers
    for(id in het) {
      p = pat(s, id)
      m = mat(s, id)
      diff = p != m
      eq1 = eq1 & diff & (a1 == p | a1 == m)
      eq2 = eq2 & diff & (a2 == p | a2 == m)
    }
    
    # Remaining carriers
    for(id in car) {
      p = pat(s, id)
      m = mat(s, id)
      eq1 = eq1 & (a1 == p | a1 == m)
      eq2 = eq2 & (a2 == p | a2 == m)
    }
    
    # Reduce sim matrix
    s = s[eq1 | eq2, , drop = FALSE]
    
    # Noncarriers (TODO: can probably be optimised)
    if(length(non)) {
      keep = rep(TRUE, nrow(s))
        
      # Which rows used a1 vs a2
      keep1 = eq1[eq1 | eq2]
      keep2 = eq2[eq1 | eq2]
      
      s1 = s[keep1, , drop = FALSE]
      s2 = s[keep2, , drop = FALSE]

      # shared alleles
      a1 = pat(s1, id1)
      a2 = mat(s2, id1)
      
      for(id in non) {
        p1 = pat(s1, id)
        m1 = mat(s1, id)
        p2 = pat(s2, id)
        m2 = mat(s2, id)
        
        keep[keep1] = keep[keep1] & a1 != p1 & a1 != m1
        keep[keep2] = keep[keep2] & a2 != p2 & a2 != m2
      }
      
      s = s[keep, , drop = FALSE]
    }
    
    # Merge adjacent segments
    if(merge)
      s = mergeSegments(s, checkAdjacency = TRUE)
    
    # Apply length cutoff
    if(cutoff > 0)
      s = s[s[, 'length'] > cutoff, , drop = FALSE]
    
    # Result for this sim
    s
  })
  
  if(single) 
    res = res[[1]]
  
  res
}

