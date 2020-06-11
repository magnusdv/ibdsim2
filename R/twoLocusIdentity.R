#' Estimate one-locus and two-locus identity coefficients
#'
#' Estimate by simulation the 9 (condensed) identity coefficients of a pair
#' pedigree members, and also the 9*9 matrix of _two-locus_ identity
#' coefficients at a given recombination rate. (For a description of the
#' identity states, see here: <https://github.com/magnusdv/ribd>.
#'
#' @inheritParams estimateIBD
#'
#' @return `estimateOneLocusIdentity()` returns a numeric vector of length 9,
#'   with the estimated coefficients.
#'
#'   `estimateTwoLocusIBD()` returns symmetric, numerical 9*9 matrix, with the
#'   estimated two-locus coefficients.
#'
#' @examples
#' x = fullSibMating(1)
#' rho = 0.25
#' Nsim = 100 # (increase for more accurate estimates!)
#'
#' estimateTwoLocusIdentity(x, ids = 5:6, rho = rho, Nsim = Nsim)
#'
#' @name estimateIdentity
NULL


#' @rdname estimateIdentity
#' @export
estimateTwoLocusIdentity = function(x, ids, rho = NULL, cM = NULL, Nsim, 
                                    Xchrom = FALSE, verbose = FALSE, ...) {
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
    m1 = estimateOneLocusIdentity(x, ids, Nsim = Nsim, Xchrom = Xchrom, verbose = verbose, ...)
    m2 = estimateOneLocusIdentity(x, ids, Nsim = Nsim, Xchrom = Xchrom, verbose = FALSE) # dont repeat verbose output
    res = outer(m1, m2)
    return(res)
  }
  
  # Define map
  map = uniformMap(cM = cM, chrom = if (Xchrom) "X" else 1)
  
  # Simulate data
  simdata = ibdsim(x, N = Nsim, ids = ids, map = map, model = "haldane", verbose = verbose, ...)
  
  # For each sim, extract first and last rows entry of column "IBD".
  sigma.list = lapply(simdata, function(a) a[c(1, nrow(a)), 'Sigma'])
  
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
  
  # If X & male(s)
  if (Xchrom) {
    sex = getSex(x, ids)
    if(sex[1] == 1 && sex[2] == 1)
      res[3:9, ] = res[, 3:9] = NA
    else if(sex[1] == 1 && sex[2] == 2)
      res[5:9, ] = res[, 5:9] = NA
    else if(sex[1] == 2 && sex[2] == 1)
      res[c(3:4, 7:9), ] = res[, c(3:4, 7:9)] = NA  
  }
  
  if (verbose) 
    cat("Total time used:", (proc.time() - st)[["elapsed"]], "seconds.\n")
  
  res / Nsim
}

#' @rdname estimateIdentity
#' @export
estimateOneLocusIdentity = function(x, ids, Nsim, Xchrom = FALSE, verbose = FALSE, ...) {
  st = proc.time()
  
  if (anyNA(match(ids, labels(x))))
    stop2("Unknown ID label: ", setdiff(ids, labels(x)))
  
  # Define map of length 0
  map = uniformMap(cM = 0, chrom = if (Xchrom) "X" else 1)
  
  # Simulate data
  simdata = ibdsim(x, N = Nsim, ids = ids, map = map, model = "haldane", verbose = verbose, ...)
  
  # For each sim, extract first entry of column "Sigma".
  ibdres = vapply(simdata, function(a) a[1, 'Sigma'], FUN.VALUE = 1)
  
  # Frequency table
  res = table(factor(ibdres, levels = 1:9, labels = paste0("state", 1:9)))
  
  # table -> vector
  res = structure(as.vector(res), names = names(res))
  
  # If X & male(s)
  if (Xchrom) {
    sex = getSex(x, ids)
    if(sex[1] == 1 && sex[2] == 1)
      res[3:9] = NA
    else if(sex[1] == 1 && sex[2] == 2)
      res[5:9] = NA
    else if(sex[1] == 2 && sex[2] == 1)
      res[c(3:4, 7:9)] = NA  
  }
  
  if (verbose) 
    cat("Total time used:", (proc.time() - st)[["elapsed"]], "seconds.\n")
  
  res / Nsim
}
