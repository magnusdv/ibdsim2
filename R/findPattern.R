#' Find specific IBD patterns
#'
#' Find segments satisfying a particular pattern of IBD sharing, in a list of
#' IBD simulations
#'
#' For each simulation, this function extracts the subset of rows satisfying the
#' allele sharing specifyed by `pattern`. That is, segments where some allele
#' occur in all of `pattern$carriers`, but in none of `pattern$noncarriers`.
#'
#' @param sims A `genomeSim` object, or a list of such. Typically made by
#'   [ibdsim()].
#' @param pattern A list of two vectors names `carriers` and `noncarriers`.
#'
#' @return A matrix (if `sims` is a single `genomeSim` object), 7or a list of
#'   matrices.
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
#' # Visual confirmation:
#' haploDraw(x, s1, margin = c(5,3,3,3))
#' 
#' @export
findPattern = function(sims, pattern) {
  
  if(single <- is.matrix(sims))
    sims = list(sims)
  
  allids = c(pattern$carriers, pattern$noncarriers)
  
  # id columns of carriers
  ca = lapply(pattern$carriers, paste0, c(":p", ":m"))
  nc = lapply(pattern$noncarriers, paste0, c(":p", ":m"))
  
  if(length(ca) == 0) 
    return(sims)
  
  # Columns names of first carrier
  ca1 = ca[[1]]
  
  res = lapply(sims, function(s) {
    s = segmentSummary(s, allids)
    eq1 = eq2 = rep(T, length = nrow(s))
    
    # The possible alleles to be shared
    a1 = s[, ca1[1]]
    a2 = s[, ca1[2]]
    
    # Find rows satisfying carriers
    for(idcols in ca[-1]) {
      patcol = s[, idcols[1]]
      matcol = s[, idcols[2]]
      eq1 = eq1 & (a1 == patcol | a1 == matcol)
      eq2 = eq2 & (a2 == patcol | a2 == matcol)
    }
    
    s = s[eq1 | eq2, , drop = FALSE]
    
    if(length(nc) == 0)
      return(s)
    
    ### Test noncarriers (reduce first)
    keep = rep(T, nrow(s))
    
    # which rows used a1 vs a2
    keep1 = eq1[eq1 | eq2]
    keep2 = eq2[eq1 | eq2]
    
    # shared alleles
    a1 = s[keep1, ca1[1]]
    a2 = s[keep2, ca1[2]]
    
    for(idcols in nc) {
      pat1 = s[keep1, idcols[1]]
      mat1 = s[keep1, idcols[2]]
      pat2 = s[keep2, idcols[1]]
      mat2 = s[keep2, idcols[2]]
      
      keep[keep1] = keep[keep1] & a1 != pat1 & a1 != mat1
      keep[keep2] = keep[keep2] & a2 != pat2 & a2 != mat2
    }
    
    s[keep, , drop = FALSE]
  })
  
  if(single) 
    res = res[[1]]
  
  res
}