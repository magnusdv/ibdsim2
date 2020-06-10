#' Estimate one-locus and two-locus inbreeding coefficients
#'
#' Estimate by simulation the inbreeding coefficient of any
#' pedigree member, and also the _two-locus_ inbreeding
#' coefficient at a given recombination rate.
#'
#' @inheritParams estimateIBD
#' @param id A single ID label
#' 
#' @return A single numeric
#'
#' @examples
#' x = halfCousinPed(0, child = TRUE)
#' rho = 0.25
#' Nsim = 100 # Increase!
#'
#' estimateTwoLocusInbreeding(x, id = 6, rho = rho, Nsim = Nsim)
#'
#' @name estimateInbreeding
NULL


#' @rdname estimateInbreeding
#' @export
estimateTwoLocusInbreeding = function(x, id, rho = NULL, cM = NULL, Nsim, 
                                      Xchrom = FALSE, verbose = FALSE, ...) {
  st = proc.time()
  
  if (!id %in%  labels(x))
    stop2("Unknown ID label: ", id)
  
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
    m1 = estimateOneLocusInbreeding(x, id, Nsim = Nsim, Xchrom = Xchrom, verbose = verbose, ...)
    m2 = estimateOneLocusInbreeding(x, id, Nsim = Nsim, Xchrom = Xchrom, verbose = FALSE) # dont repeat verbose output
    res = m1 * m2
    return(res)
  }
  
  # Define map
  map = uniformMap(cM = cM, chromosome = if (Xchrom) 23 else 1)
  
  # Simulate data
  simdata = ibdsim(x, sims = Nsim, ids = id, map = map, model = "haldane", verbose = verbose, ...)
  
  # For each sim, extract first and last rows entry of column "IBD".
  f2 = lapply(simdata, function(a) {
    a[1, 5] == a[1, 6] && a[nrow(a), 5] == a[nrow(a), 6] 
  })
  
  # Shape list of observations into wide matrix
  f2 = unlist(f2)
  
  if (verbose) 
    cat("Total time used:", (proc.time() - st)[["elapsed"]], "seconds.\n")
  
  mean(f2)
}

#' @rdname estimateInbreeding
#' @export
estimateOneLocusInbreeding = function(x, id, Nsim, Xchrom = FALSE, verbose = FALSE, ...) {
  st = proc.time()
  
  if (!id %in%  labels(x))
    stop2("Unknown ID label: ", id)
  
  # Define map of length 0
  map = uniformMap(cM = 0, chromosome = if (Xchrom) 23 else 1)
  
  # Simulate data
  simdata = ibdsim(x, sims = Nsim, ids = id, map = map, model = "haldane", verbose = verbose, ...)
  
  # For each sim, check for autozygosity
  f2 = vapply(simdata, function(a) {a[1, 5] == a[1, 6]}, FUN.VALUE = TRUE)
  
  if (verbose) 
    cat("Total time used:", (proc.time() - st)[["elapsed"]], "seconds.\n")
  
  mean(f2)
}

