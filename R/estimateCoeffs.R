#' Estimation of one- and two-locus relatedness coefficients
#'
#' Estimate by simulation various relatedness coefficients, and two-locus
#' versions of the same coefficients, for a given recombination rate. The
#' current implementation covers inbreeding coefficients, kinship coefficients,
#' IBD (kappa) coefficients between noninbred individuals, and condensed
#' identity coefficients. These functions are primarily meant as tools for
#' validating exact algorithms, e.g., as implemented in the `ribd` package.
#'
#' In the following, let L1 and L2 denote two arbitrary autosomal loci with
#' recombination rate \eqn{\rho}, and let A and B be members of the pedigree
#' `x`.
#'
#' The *two-locus inbreeding coefficient* \eqn{f_2(\rho)} of A is defined as the
#' probability that A is autozygous at both L1 and L2 simultaneously.
#'
#' The *two-locus kinship coefficient* \eqn{\phi_2(\rho)} of A and B is defined
#' as the probability that a random gamete emitted from A, and a random gamete
#' emitted from B, contain IBD alleles at both L1 and L2.
#'
#' The *two-locus kappa coefficient* \eqn{\kappa_{ij}(\rho)}, for \eqn{i,j =
#' 0,1,2}, of noninbred A and B, is the probability that A and B share exactly
#' `i` alleles IBD at L1, and exactly `j` alleles IBD at L2.
#'
#' The *two-locus identity coefficient* \eqn{\Delta_{ij}}, \eqn{i,j = 1,...,9}
#' is defined for any (possibly inbred) A and B, as the probability that A and B
#' are in identity state `i` at L1, and state `j` at L2. This uses the
#' conventional ordering of the nine condensed identity states. For details, see
#' for instance the [GitHub page of the `ribd`
#' package](https://github.com/magnusdv/ribd).
#'
#' @param x A pedigree in the form of a [pedtools::ped()] object.
#' @param id,ids A vector of one or two ID labels.
#' @param rho A scalar in the interval `[0, 0.5]`: the recombination fraction
#'   between the two loci, converted to centiMorgans using Haldane's map function:
#'   cM = -50 * log(1 - 2 * rho). Either `rho` or `cM` (but not both) must be
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
#' @return `estimateInbreeding()`: a single probability.
#'
#'   `estimateTwoLocusInbreeding()`: a single probability.
#'
#'   `estimateKappa()`: a numeric vector of length 3, with the estimated
#'   \eqn{\kappa} coefficients.
#'
#'   `estimateTwoLocusKappa()`: a symmetric, numerical 3*3 matrix, with the
#'   estimated values of \eqn{\kappa_{ij}}, for \eqn{i,j = 0,1,2}.
#'
#'   `estimateIdentity()`: a numeric vector of length 9, with the estimated
#'   identity coefficients.
#'
#'   `estimateTwoLocusIdentity()`: a symmetric, numerical 9*9 matrix, with the
#'   estimated values of \eqn{\Delta_{ij}}, for \eqn{i,j = 1,...,9}.
#'
#'
#' @examples
#'
#' ############################
#' ### Two-locus inbreeding ###
#' ############################
#'
#' x = cousinPed(0, child = TRUE)
#' rho = 0.25
#' Nsim = 10 # Increase!
#' estimateTwoLocusInbreeding(x, id = 5, rho = rho, Nsim = Nsim, seed = 123)
#'
#' ########################################
#' ### Two-locus kappa:                 ###
#' ### Grandparent vs half sib vs uncle ###
#' ########################################
#'
#' # These are indistinguishable with unlinked loci, see e.g.
#' # pages 182-183 in Egeland, Kling and Mostad (2016).
#' # In the following, each simulation approximation is followed
#' # by its exact counterpart.
#'
#' rho = 0.25; R = .5 * (rho^2 + (1-rho)^2)
#' Nsim = 10 # Should be increased to at least 10000
#'
#' # Grandparent/grandchild
#' G = linearPed(2); G.ids = c(1,5); # plot(G, hatched = G.ids)
#' estimateTwoLocusKappa(G, G.ids, rho = rho, Nsim = Nsim, seed = 123)[2,2]
#' .5*(1-rho) # exact
#'
#' # Half sibs
#' H = halfSibPed(); H.ids = c(4,5); # plot(H, hatched = H.ids)
#' estimateTwoLocusKappa(H, H.ids, rho = rho, Nsim = Nsim, seed = 123)[2,2]
#' R # exact
#'
#' # Uncle
#' U = cousinPed(0, removal = 1); U.ids = c(3,6); # plot(U, hatched = U.ids)
#' estimateTwoLocusKappa(U, U.ids, rho = rho, Nsim = Nsim, seed = 123)[2,2]
#' (1-rho) * R + rho/4 # exact
#'
#' # Exact calculations by ribd:
#' # ribd::twoLocusIBD(G, G.ids, rho = rho, coefs = "k11")
#' # ribd::twoLocusIBD(H, H.ids, rho = rho, coefs = "k11")
#' # ribd::twoLocusIBD(U, U.ids, rho = rho, coefs = "k11")
#'
#' ##########################
#' ### Two-locus Jacquard ###
#' ##########################
#'
#' x = fullSibMating(1)
#' rho = 0.25
#' Nsim = 10 # (increase to at least 10000)
#'
#' estimateTwoLocusIdentity(x, ids = 5:6, rho = rho, Nsim = Nsim, seed = 123)
#'
#' # Exact by ribd:
#' # ribd::twoLocusIdentity(x, ids = 5:6, rho = rho)
#'
#' @name estimateCoeffs
NULL


#' @rdname estimateCoeffs
#' @export
estimateInbreeding = function(x, id, Nsim, Xchrom = FALSE, verbose = FALSE, ...) {
  st = proc.time()
  
  # Define map of length 0
  map = uniformMap(cM = 0, chrom = if (Xchrom) "X" else 1)
  
  # Simulate data
  simdata = ibdsim(x, N = Nsim, ids = id, map = map, model = "haldane", verbose = verbose, ...)
  
  # For each sim, extract autozygosity status
  f2 = vapply(simdata, function(a) a[[1, "Aut"]], FUN.VALUE = 1)
  
  if (verbose) 
    cat("Total time used:", (proc.time() - st)[["elapsed"]], "seconds.\n")
  
  mean(f2)
}

#' @rdname estimateCoeffs
#' @export
estimateTwoLocusInbreeding = function(x, id, rho = NULL, cM = NULL, Nsim, 
                                      Xchrom = FALSE, verbose = FALSE, ...) {
  st = proc.time()
  
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
    m1 = estimateInbreeding(x, id, Nsim = Nsim, Xchrom = Xchrom, verbose = verbose, ...)
    m2 = estimateInbreeding(x, id, Nsim = Nsim, Xchrom = Xchrom, verbose = FALSE) # dont repeat verbose output
    res = m1 * m2
    return(res)
  }
  
  # Define map
  map = uniformMap(cM = cM, chrom = if (Xchrom) "X" else 1)
  
  # Simulate data
  simdata = ibdsim(x, N = Nsim, ids = id, map = map, model = "haldane", verbose = verbose, ...)
  
  # For each sim, check if autozygous at both ends
  f2 = vapply(simdata, function(a) a[[1, "Aut"]] == 1 && a[[nrow(a), "Aut"]] == 1, FUN.VALUE = FALSE)
  
  if (verbose) 
    cat("Total time used:", (proc.time() - st)[["elapsed"]], "seconds.\n")
  
  mean(f2)
}


#' @rdname estimateCoeffs
#' @export
estimateKinship = function(x, ids, Nsim, Xchrom = FALSE, verbose = FALSE, ...) {
  st = proc.time()
  
  # Define map of length 0
  map = uniformMap(cM = 0, chrom = if (Xchrom) "X" else 1)
  
  # Simulate data
  simdata = ibdsim(x, N = Nsim, ids = ids, map = map, model = "haldane", verbose = verbose, ...)
  
  # For each sim, extract entry in column "Sigma".
  jacq = vapply(simdata, function(a) a[[1, 'Sigma']], FUN.VALUE = 1)
  
  if (verbose) 
    cat("Total time used:", (proc.time() - st)[["elapsed"]], "seconds.\n")
  
  # Coefficients in the relation phi = weights * jacq
  wei = c(1, 0, .5, 0, .5, 0, .5, .25, 0)
  
  mean(wei[jacq])
}


#' @rdname estimateCoeffs
#' @export
estimateTwoLocusKinship = function(x, ids, rho = NULL, cM = NULL, Nsim, 
                                 Xchrom = FALSE, verbose = FALSE, ...) {
  stop2("Not implemented yet; please contact the developer if you need this")
  st = proc.time()
  
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
    m1 = estimateKappa(x, ids, Nsim = Nsim, Xchrom = Xchrom, verbose = verbose, ...)
    m2 = estimateKappa(x, ids, Nsim = Nsim, Xchrom = Xchrom, verbose = FALSE, ...) # don't repeat verbose output
    res = outer(m1, m2)
    return(res)
  }
  
  # Define map
  map = uniformMap(cM = cM, chrom = if (Xchrom) "X" else 1)
  
  # Simulate data
  simdata = ibdsim(x, N = Nsim, ids = ids, map = map, model = "haldane", verbose = verbose, ...)
  
  # For each sim, extract Sigma states at both ends
  jacq = lapply(simdata, function(a) a[c(1, nrow(a)), "Sigma"])
  
  # TODO: Finish this. Needs coefficient matrix expressed by rho, 1-rho, R = rho^2 + (1-rho)^2 
}


#' @rdname estimateCoeffs
#' @export
estimateKappa = function(x, ids, Nsim, Xchrom = FALSE, verbose = FALSE, ...) {
  st = proc.time()
  
  # Define map of length 0
  map = uniformMap(cM = 0, chrom = if (Xchrom) "X" else 1)
  
  # Simulate data
  simdata = ibdsim(x, N = Nsim, ids = ids, map = map, model = "haldane", verbose = verbose, ...)
  
  # For each sim, extract entry in column "IBD".
  ibdres = vapply(simdata, function(a) a[[1, 'IBD']], FUN.VALUE = 1)
  
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


#' @rdname estimateCoeffs
#' @export
estimateTwoLocusKappa = function(x, ids, rho = NULL, cM = NULL, Nsim, 
                                 Xchrom = FALSE, verbose = FALSE, ...) {
  st = proc.time()

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
    m1 = estimateKappa(x, ids, Nsim = Nsim, Xchrom = Xchrom, verbose = verbose, ...)
    m2 = estimateKappa(x, ids, Nsim = Nsim, Xchrom = Xchrom, verbose = FALSE, ...) # don't repeat verbose output
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


#' @rdname estimateCoeffs
#' @export
estimateIdentity = function(x, ids, Nsim, Xchrom = FALSE, verbose = FALSE, ...) {
  st = proc.time()
  
  # Define map of length 0
  map = uniformMap(cM = 0, chrom = if (Xchrom) "X" else 1)
  
  # Simulate data
  simdata = ibdsim(x, N = Nsim, ids = ids, map = map, model = "haldane", verbose = verbose, ...)
  
  # For each sim, extract entry in column "Sigma".
  ibdres = vapply(simdata, function(a) a[[1, 'Sigma']], FUN.VALUE = 1)
  
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


#' @rdname estimateCoeffs
#' @export
estimateTwoLocusIdentity = function(x, ids, rho = NULL, cM = NULL, Nsim, 
                                    Xchrom = FALSE, verbose = FALSE, ...) {
  st = proc.time()
  
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
    m1 = estimateIdentity(x, ids, Nsim = Nsim, Xchrom = Xchrom, verbose = verbose, ...)
    m2 = estimateIdentity(x, ids, Nsim = Nsim, Xchrom = Xchrom, verbose = FALSE) # dont repeat verbose output
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