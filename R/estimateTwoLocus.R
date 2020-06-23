#' Estimation of one- and two-locus relatedness coefficients
#'
#' Estimate by simulation various relatedness coefficients, and two-locus
#' versions of the same coefficients, for a given recombination rate. The
#' current implementation covers inbreeding coefficients, pairwise IBD (kappa)
#' coefficients and pairwise condensed identity coefficients. These functions
#' are primarily meant as tools for validating exact algorithms, e.g., as
#' implemented in the `ribd` package.
#'
#' In the following, let L1 and L2 denote two arbitrary autosomal loci with
#' recombination rate \eqn{\rho}.
#'
#' The *two-locus inbreeding coefficient* \eqn{f_2(\rho)} of a pedigree member
#' is defined as the probability of being autozygous at both L1 and L2
#' simultaneously.
#'
#' The *two-locus IBD coefficients* \eqn{\kappa_{ij}(\rho)}, for \eqn{i,j =
#' 0,1,2}, of individuals A and B are defined by as the probability that A and B
#' share i alleles IBD at L1, and j alleles IBD at L2.
#'
#' The *two-locus identity coefficients* \eqn{\Delta_{ij}}, \eqn{i,j = 1,...,9}
#' are defined similarly to the two-locus kappa above. For a description of the
#' identity states, see e.g., <https://github.com/magnusdv/ribd>.
#'
#' @param x A pedigree in the form of a [pedtools::ped()] object.
#' @param id,ids A vector of one or two ID labels.
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
#' @return `estimateOneLocusInbreeding()`: a single probability.
#'
#'   `estimateTwoLocusInbreeding()`: a single probability.
#'
#'   `estimateOneLocusIBD()`: a numeric vector of length 3, with the estimated
#'   \eqn{\kappa} coefficients.
#'
#'   `estimateTwoLocusIBD()`: a symmetric, numerical 3*3 matrix, with the
#'   estimated values of \eqn{\kappa_{ij}}, for \eqn{i,j = 0,1,2}.
#'
#'   `estimateOneLocusIdentity()`: a numeric vector of length 9, with the
#'   estimated identity coefficients.
#'
#'   `estimateTwoLocusIBD()`: symmetric, numerical 9*9 matrix, with the
#'   estimated two-locus identity coefficients.
#'
#'
#' @examples
#'
#' ### Two-locus inbreeding ###
#'
#' x = cousinPed(0, child = TRUE)
#' rho = 0.25
#' Nsim = 10 # Increase!
#' estimateTwoLocusInbreeding(x, id = 5, rho = rho, Nsim = Nsim)
#'
#'
#' ### Two-locus IBD: Grandparent vs half sib vs uncle ###
#'
#' # These are indistinguishable with unlinked loci, see e.g.
#' # pages 182-183 in Egeland, Kling and Mostad (2016).
#' # Each simulations followed by exact counterpart.
#'
#' rho = 0.25; R = .5 * (rho^2 + (1-rho)^2)
#' Nsim = 10 # Should be increased to at least 10000
#'
#' # Grandparent/grandchild
#' G = linearPed(2); G.ids = c(1,5); #plot(G, shaded = G.ids)
#' estimateTwoLocusIBD(G, G.ids, rho = rho, Nsim = Nsim)[2,2]
#' .5*(1-rho) # exact
#'
#' # Half sibs
#' H = halfSibPed(); H.ids = c(4,5); # plot(H, shaded = H.ids)
#' estimateTwoLocusIBD(H, H.ids, rho = rho, Nsim = Nsim)[2,2]
#' R # exact
#'
#' # Uncle
#' U = cousinPed(0, removal = 1); U.ids = c(3,6); # plot(U, shaded = U.ids)
#' estimateTwoLocusIBD(U, U.ids, rho = rho, Nsim = Nsim)[2,2]
#' (1-rho) * R + rho/4 # exact
#'
#'
#' ### Two-locus Jacquard ###
#' x = fullSibMating(1)
#' rho = 0.25
#' Nsim = 10 # (increase for more accurate estimates!)
#'
#' estimateTwoLocusIdentity(x, ids = 5:6, rho = rho, Nsim = Nsim)
#'
#' @name estimateTwoLocus
NULL


#' @rdname estimateTwoLocus
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
  map = uniformMap(cM = cM, chrom = if (Xchrom) "X" else 1)
  
  # Simulate data
  simdata = ibdsim(x, N = Nsim, ids = id, map = map, model = "haldane", verbose = verbose, ...)
  
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

#' @rdname estimateTwoLocus
#' @export
estimateOneLocusInbreeding = function(x, id, Nsim, Xchrom = FALSE, verbose = FALSE, ...) {
  st = proc.time()
  
  if (!id %in%  labels(x))
    stop2("Unknown ID label: ", id)
  
  # Define map of length 0
  map = uniformMap(cM = 0, chrom = if (Xchrom) "X" else 1)
  
  # Simulate data
  simdata = ibdsim(x, N = Nsim, ids = id, map = map, model = "haldane", verbose = verbose, ...)
  
  # For each sim, check for autozygosity
  f2 = vapply(simdata, function(a) {a[1, 5] == a[1, 6]}, FUN.VALUE = TRUE)
  
  if (verbose) 
    cat("Total time used:", (proc.time() - st)[["elapsed"]], "seconds.\n")
  
  mean(f2)
}


#' @rdname estimateTwoLocus
#' @export
estimateTwoLocusIBD = function(x, ids, rho = NULL, cM = NULL, Nsim, 
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
    m1 = estimateOneLocusIBD(x, ids, Nsim = Nsim, Xchrom = Xchrom, verbose = verbose, ...)
    m2 = estimateOneLocusIBD(x, ids, Nsim = Nsim, Xchrom = Xchrom, verbose = FALSE, ...) # dont repeat verbose output
    res = outer(m1, m2)
    return(res)
  }

  # Define map
  map = uniformMap(cM = cM, chrom = if (Xchrom) "X" else 1)
  
  # Simulate data
  simdata = ibdsim(x, N = Nsim, ids = ids, map = map, model = "haldane", verbose = verbose, ...)

  # For each sim, extract first and last rows entry of column "IBD".
  ibd.list = lapply(simdata, function(a) a[c(1, nrow(a)), 'IBD'])
  
  # Shape IBD list into matrix
  ibd.mat = unlist(ibd.list)
  dim(ibd.mat) = c(2, Nsim)
  
  # Frequency table
  labs = paste0("ibd", 0:2)
  res = table(factor(ibd.mat[1, ], levels = 0:2, labels = labs), 
              factor(ibd.mat[2, ], levels = 0:2, labels = labs))
  
  # table -> matrix
  res = unclass(res) 
  names(dimnames(res)) = NULL
  
  # If X, set IBD2 = NA for males
  if (Xchrom && any(getSex(x, ids) == 1))
    res['ibd2', ] = res[, 'ibd2'] = NA
  
  if (verbose) 
    cat("Total time used:", (proc.time() - st)[["elapsed"]], "seconds.\n")
  
  res / Nsim
}

#' @rdname estimateTwoLocus
#' @export
estimateOneLocusIBD = function(x, ids, Nsim, Xchrom = FALSE, verbose = FALSE, ...) {
  st = proc.time()
  
  if (anyNA(match(ids, labels(x))))
    stop2("Unknown ID label: ", setdiff(ids, labels(x)))
  
  # Define map of length 0
  map = uniformMap(cM = 0, chrom = if (Xchrom) "X" else 1)
  
  # Simulate data
  simdata = ibdsim(x, N = Nsim, ids = ids, map = map, model = "haldane", verbose = verbose, ...)
  
  # For each sim, extract first and last rows entry of column "IBD".
  ibdres = vapply(simdata, function(a) a[, 'IBD'], FUN.VALUE = 1)
  
  # Frequency table
  res = table(factor(ibdres, levels = 0:2, labels = paste0("ibd", 0:2)))
  
  # table -> vector
  res = structure(as.vector(res), names = names(res))
  
  # If X & male(s), set IBD2 = NA
  if (Xchrom && any(getSex(x, ids) == 1)) 
    res['ibd2'] = NA
  
  if (verbose) 
    cat("Total time used:", (proc.time() - st)[["elapsed"]], "seconds.\n")
  
  res / Nsim
}


#' @rdname estimateTwoLocus
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

#' @rdname estimateTwoLocus
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


