#' Probability of zero IBD
#'
#' Estimate the probability of no IBD sharing in a pairwise relationship.
#'
#' @param sims A list of genome simulations, as output by [ibdsim()].
#' @param ids A vector with two ID labels. If NULL (default), these are deduced
#'   from the `sims` object.
#' @param threshold A nonnegative number (default:0). Only IBD segments longer
#'   than this are included in the computation.
#' @param unit The unit of measurement for `threshold`: Either "mb" or "cm"
#'   (default).
#'
#' @return A list with the following two entries:
#'
#'   * `zeroprob`: The fraction of `sims` in which `ids` have no IBD sharing
#'
#'   * `stErr`: The standard error of `zeroprob`
#'
#' @examples
#' ###
#' # The following example computes the probability of
#' # no IBD sharing between a pair of fourth cousins.
#' # We also show how the probability is affected by
#' # truncation, i.e., ignoring short segments.
#' ###
#'
#' # Define the pedigree
#' x = cousinPed(4)
#' cous = leaves(x)
#'
#' # Simulate (increase N!)
#' s = ibdsim(x, N = 10)
#'
#' # Probability of zero ibd segments. (By default all segs are used)
#' zeroIBD(s, ids = cous)
#'
#' # Re-compute with nonzero threshold
#' zeroIBD(s, ids = cous, threshold = 1, unit = "cm")
#' zeroIBD(s, ids = cous, threshold = 1, unit = "mb")
#'
#' @export
zeroIBD = function(sims, ids = NULL, threshold = 0, unit = "cm") {
  if(!is.numeric(threshold) || length(threshold) != 1 || threshold < 0)
    stop2("`threshold` must be a nonnegative number")
  
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
  
  # Names of start/end columns
  startCol = switch(unit, mb = "startMB", cm = "startCM")
  endCol = switch(unit, mb = "endMB", cm = "endCM")
  
  # Summarise each simulation
  ibdCount = lapply(sims, function(s) {
    
    if(length(idsims) > 2 || !"IBD" %in% colnames(s)) {
      s0 = alleleFlow(s, ids = ids, addState = TRUE)
      s = mergeSegments(s0, by = "IBD")
    }
    
    len = s[, endCol] - s[, startCol]
    ibdstate = s[, 'IBD']
    
    sum(ibdstate > 0 & len >= threshold)
  })
  
  # Fraction (and standard error) of sims with 0 segments
  zeroprob = mean(unlist(ibdCount) == 0)
  stErr = sqrt(zeroprob*(1 - zeroprob)/length(ibdCount))
  
  list(zeroprob = zeroprob, stErr = stErr)
}
