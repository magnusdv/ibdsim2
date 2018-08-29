#' Probability of zero ibd
#'
#' Estimate the probability of no ibd sharing in a pairwise relationship.
#'
#' @param sim A list of genome simulations, as output by [ibdsim()].
#' @param id.pair A vector of length 2, with ID labels of the two individuals in
#'   question.
#' @param truncate A vector of positive real numbers. Only IBD segments longer than this are
#'   included in the computation. If `truncate` has more than one
#'   element, a separate estimate is provided for each value. The default (`truncate = 0`) is to
#'   include all segments.
#'   
#' @return A data.frame with tree numeric columns:
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
#' x = pedtools::cousinPed(4)
#' cousins = pedtools::leaves(x)
#' 
#' # Simulate (increase 'sims'!)
#' s = ibdsim(x, sims = 7)
#' 
#' # Probability of zero ibd segments. (By default all segs are used)
#' zero_ibd(s, id.pair = cousins)
#' 
#' # Re-compute the probability with several truncation levels
#' truncate = 0:50
#' zp = zero_ibd(s, id.pair = cousins, truncate = truncate)
#' 
#' plot(truncate, zp$zeroprob)
#' 
#' @export
zero_ibd = function(sim, id.pair, truncate=0) {
  if(length(id.pair) != 2)
    stop2("`id.pair` must be a vector of length 2")
  if(!is.numeric(truncate) || length(truncate) == 0 || any(truncate < 0))
    stop2("`truncate` must be vector of positive numbers"))
  
  ibd_count = vapply(sim, function(s) {
    a = alleleSummary(s, ids=id.pair, ibd.status=T)
    ibdstatus = a[, 'ibd']
    len = a[, 'length']
    
    # Count IBD segments (ibd = 1 or 2) longer than each "truncate"
    segs = vapply(truncate, function(trunc) sum(ibdstatus > 0 & len >= trunc), FUN.VALUE=1)
    
  }, FUN.VALUE=numeric(length(truncate)))
  
  # Fix vapply output inconsistency.
  if(length(truncate) == 1) dim(ibd_count) = c(1, length(ibd_count))
  
  # Fraction (and standard error) of sims with 0 segments
  zeroprob = rowMeans(ibd_count == 0)
  SE = sqrt(zeroprob*(1 - zeroprob)/length(sim))
  
  data.frame(truncate=truncate, zeroprob = zeroprob, SE = SE)
}
