#' Estimate one-locus and two-locus identity coefficients
#'
#' Estimate by simulation the 9 (condensed) identity coefficients of a pair
#' pedigree members, and also the 9*9 matrix of _two-locus_ identity
#' coefficients at a given recombination rate.
#'
#' @param x A pedigree in the form of a [pedtools::ped()] object.
#' @param ids A character (or coercible to character) of length 2.
#' @param rho A scalar in the interval `[0, 0.5]`: the recombination fraction
#'   between the two loci, converted to centiMorgan using Haldanes map function:
#'   cM = -50 * log(1 - 2*rho). Either `rho` or `cM` (but not both) must be
#'   non-NULL.
#' @param cM A non-negative number: the genetic distance between the two loci,
#'   given in centiMorgans. Either `rho` or `cM` (but not both) must be
#'   non-NULL.
#' @param Nsim The number of simulations.
#' @param Xchrom A logical indicating if the loci are X-linked (if TRUE) or
#'   autosomal (FALSE).
#' @param verbose A logical.
#' @param ... Further arguments passed on to [ibdsim()], e.g. `seed`.
#'
#' @return `estimateOneLocusIdentity()` returns a numeric vector of length 9, with
#'   the estimated coefficients.
#'
#'   `estimateTwoLocusIBD()` returns symmetric, numerical 9*9 matrix, with the
#'   estimated two-locus coefficients.
#'
#' @examples
#' x = fullSibMating(1)
#' rho = 0.25
#' Nsim = 100 # Increase!
#'
#' estimateTwoLocusIdentity(x, ids = 5:6, rho = rho, Nsim = Nsim)
#'
#' @name estimateIdentity
NULL

#' @rdname estimateIdentity
#' @export
estimateTwoLocusIdentity = function(x, ids, rho = NULL, cM = NULL, Nsim, 
                                    Xchrom = F, verbose = F, ...) {
  st = proc.time()
  
  if (anyNA(match(ids, labels(x))))
    stop2("Unknown ID label: ", setdiff(ids, labels(x)))
  
  if (is.null(cM) + is.null(rho) != 1) 
    stop2("Exactly one of the parameters `cM` and `rho` must be non-NULL")
  
  if (!is.null(rho)) {
    if (rho < 0 || rho > 0.5) 
      stop2("Recombination rate `rho` is outside of interval [0, 0.5]: ", rho)
    cM = -50 * log(1 - 2 * rho)
  }
  if (verbose) cat("Locus distance:", cM, "centiMorgan\n")
  
  if (cM == Inf) {
    if (verbose) cat("Analysing unlinked loci.\n")
    m1 = estimateOneLocusIdentity(x, ids, Nsim = Nsim, verbose = verbose, ...)
    m2 = estimateOneLocusIdentity(x, ids, Nsim = Nsim, ...) # dont repeat verbose output
    res = outer(m1, m2)
    return(res)
  }
  
  # Define map
  map = uniformMap(cM = cM, chromosome = if (Xchrom) 23 else 1)
  
  # Simulate data
  simdata = ibdsim(x, map = map, sims = Nsim, model = "haldane", verbose = verbose, ...)
  
  # For each sim, extract first and last rows entry of column "IBD".
  sumFUN = if (Xchrom) alleleSummaryX else alleleSummary
  sigma.list = lapply(simdata, function(s) {
    a = sumFUN(s, ids)
    a[c(1, nrow(a)), 'Sigma']
  })
  
  # Shape list of observations into wide matrix
  sigma.mat = unlist(sigma.list)
  dim(sigma.mat) = c(2, Nsim)
  
  # Frequency table
  labs = paste0("state", 1:9)
  res = table(factor(sigma.mat[1, ], levels = 1:9, labels = labs), 
              factor(sigma.mat[2, ], levels = 1:9, labels = labs))
  
  # table -> matrix
  res = unclass(res) 
  names(dimnames(res)) = NULL
  
  if (verbose) 
    cat("Total time used:", (proc.time() - st)[["elapsed"]], "seconds.\n")
  
  res / Nsim
}

#' @rdname estimateIdentity
#' @export
estimateOneLocusIdentity = function(x, ids, Nsim, Xchrom = F, verbose = F, ...) {
  st = proc.time()
  
  if (anyNA(match(ids, labels(x))))
    stop2("Unknown ID label: ", setdiff(ids, labels(x)))
  
  # Define map of length 0
  map = uniformMap(cM = 0, chromosome = if (Xchrom) 23 else 1)
  
  # Simulate data
  simdata = ibdsim(x, map = map, sims = Nsim, model = "haldane", verbose = verbose, ...)
  
  # For each sim, extract first entry of column "Sigma".
  sumFUN = if (Xchrom) alleleSummaryX else alleleSummary
  ibdres = vapply(simdata, function(s) sumFUN(s, ids)[1, 'Sigma'], 1)
  
  # Frequency table
  res = table(factor(ibdres, levels = 1:9, labels = paste0("state", 1:9)))
  
  # table -> vector
  res = structure(as.vector(res), names = names(res))
  
  # If X & male(s)
  # TODO!
  
  if (verbose) 
    cat("Total time used:", (proc.time() - st)[["elapsed"]], "seconds.\n")
  
  res / Nsim
}
