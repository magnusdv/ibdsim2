#' Probability of zero IBD
#'
#' Estimate the probability of no IBD sharing in a pairwise relationship.
#'
#' @param sims A list of genome simulations, as output by [ibdsim()].
#' @param ids A vector with two ID labels.
#' @param truncate A vector of positive real numbers. Only IBD segments longer than this are
#'   included in the computation. If `truncate` has more than one
#'   element, a separate estimate is provided for each value. The default (`truncate = 0`) is to
#'   include all segments.
#'   
#' @return A data frame with tree numeric columns:
#' \describe{
#'   \item{truncate}{Same as input.}
#'   \item{zeroprob}{The estimated probability of no IBD segments.}
#'   \item{SE}{The standard error of the estimate.}
#' }
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
#' cousins = leaves(x)
#' 
#' # Simulate (increase N!)
#' s = ibdsim(x, N = 10)
#' 
#' # Probability of zero ibd segments. (By default all segs are used)
#' zeroIBD(s, ids = cousins)
#' 
#' # Re-compute the probability with several truncation levels
#' truncate = 0:20
#' zp = zeroIBD(s, ids = cousins, truncate = truncate)
#' 
#' plot(truncate, zp$zeroprob, type = "b", ylim = c(0,1))
#' 
#' @export
zeroIBD = function(sims, ids, truncate = 0) {
  if(!is.numeric(truncate) || length(truncate) == 0 || any(truncate < 0))
    stop2("`truncate` must be vector of positive numbers")
  
  ibdCount = vapply(sims, function(s) {
    a = segmentSummary(s, ids, addState = TRUE)
    
    ibdstatus = a[, 'IBD']
    len = a[, 'length']
    
    # Count IBD segments (ibd = 1 or 2) longer than each "truncate"
    segs = vapply(truncate, function(trunc) sum(ibdstatus > 0 & len >= trunc), 1)
  }, 
  FUN.VALUE = numeric(length(truncate)))
  
  # Fix vapply output inconsistency.
  if(length(truncate) == 1) 
    dim(ibdCount) = c(1, length(ibdCount))
  
  # Fraction (and standard error) of sims with 0 segments
  zeroprob = rowMeans(ibdCount == 0)
  SE = sqrt(zeroprob*(1 - zeroprob)/length(sims))
  
  data.frame(truncate = truncate, zeroprob = zeroprob, SE = SE)
}
