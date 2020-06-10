#' Estimate one-locus and two-locus IBD coefficients
#'
#' Estimate by simulation the IBD coefficients \eqn{(\kappa[0], \kappa[1],
#' \kappa[2])} of two non-inbred pedigree members, and also the matrix of _two-locus_ IBD coefficients
#' \eqn{(\kappa[i,j])} at a given recombination rate.
#'
#' Alleles segregating in a pedigree are said to be _identical by descent_ (IBD) if
#' they originate from the same ancestral allele within the pedigree.
#'
#' For two non-inbred individuals A and B, their IBD coefficients of \eqn{\kappa
#' = (\kappa[0], \kappa[1], \kappa[2])} are defined as the probabilities
#' \deqn{\kappa[i] = Pr(A and B share i alleles IBD at a given autosomal locus).}

#' Similarly, the _two-locus_ IBD coefficients of A and B are defined by
#' \deqn{\kappa[i,j] = Pr(A and B share i alleles IBD at locus 1, and j alleles
#' IBD at locus 2),} where \eqn{0 \leq i,j \leq 2}{0 <= i,j <= 2}.
#'
#' While the single-locus IBD coefficients depend only on the genealogy relating
#' the two individuals, the two-locus coefficients also depend on the genetic
#' distance between the loci. In particular, if the loci are completely linked
#' (`rho = 0`, or equivalently, `cM = 0`) the IBD matrix is diagonal with
#' \eqn{\kappa[i,i] = \kappa[i]}. If the loci are completely _unlinked_ (`rho =
#' 0.5`; `cM = Inf`) then \eqn{\kappa[i,j] = \kappa[i]*\kappa[j]}. (See
#' examples.)
#'
#' X chromosomal coefficients: If either A or B is male, the IBD status at X
#' chromosomal loci cannot be 2. This is reflected by NA entries in the output.
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
#' @return `estimateOneLocusIBD()` returns a numeric vector of length 3, with
#'   the estimated values of \eqn{(\kappa[0],\kappa[1],\kappa[2]}
#'
#'   `estimateTwoLocusIBD()` returns symmetric, numerical 3*3 matrix, with the
#'   estimated values of \eqn{\kappa[i,j]}{k_i_j}, for i,j = 0,1,2.
#'
#' @examples
#' ############################
#' ### Example 1: Full siblings
#' ############################
#' x = nuclearPed(2)
#' Nsim = 100 # Should be increased substantially
#'
#' ### 1a) One-locus kappa estimates (autosomal and X):
#' k.hat = estimateOneLocusIBD(x, ids = 3:4, 
#'                             Nsim = Nsim, seed = 123)
#' k.hat.X = estimateOneLocusIBD(x, ids = 3:4, Nsim = Nsim, 
#'                               Xchrom = TRUE, seed = 123)
#'
#' ### 1b) Two-locus IBD estimation
#' # Completely linked, autosomal
#' rho = 0
#' k2.linked = estimateTwoLocusIBD(x, ids = 3:4, rho = rho, 
#'                                 Nsim = Nsim, seed = 123)
#' stopifnot(identical(diag(k2.linked), k.hat))
#'
#' # Completely unlinked, autosomal
#' rho = 0.5
#' k2.unlinked = estimateTwoLocusIBD(x, ids = 3:4, rho = rho, 
#'                                   Nsim = Nsim, seed = 123)
#' stopifnot(identical(k2.unlinked, outer(k.hat, k.hat)))
#'
#' # Recombination rate 10%, autosomal
#' rho = 0.1
#' r1 = estimateTwoLocusIBD(x, ids = 3:4, rho = rho, 
#'                          Nsim = Nsim, seed = 17)
#' 
#' # Alternatively, specify cM directly
#' cM = -50 * log(1 - 2*rho) # Haldane's map
#' r2 = estimateTwoLocusIBD(x, ids = 3:4, cM = cM, 
#'                          Nsim = Nsim, seed = 17)
#' stopifnot(identical(r1, r2))
#'
#' ### 1c) Two-locus IBD on X
#' # Completely linked
#' rho = 0
#' k2.linked.X = estimateTwoLocusIBD(x, ids = 3:4, rho = rho, Xchrom = TRUE,
#'                                   Nsim = Nsim, seed = 123)
#' stopifnot(identical(diag(k2.linked.X), k.hat.X))
#'
#' # Completely unlinked
#' rho = 0.5
#' k2.unlinked.X = estimateTwoLocusIBD(x, ids = 3:4, rho = rho, Xchrom = TRUE, 
#'                                     Nsim = Nsim, seed = 123)
#' stopifnot(identical(k2.unlinked.X, outer(k.hat.X, k.hat.X)))
#'
#' # Recombination rate 10%, X chromosome
#' rho = 0.1
#' cM = -50 * log(1 - 2*rho)
#' r1.X = estimateTwoLocusIBD(x, ids = 3:4, rho = rho, Xchrom = TRUE,
#'                            Nsim = Nsim, seed = 123)
#' r2.X = estimateTwoLocusIBD(x, ids = 3:4, cM = cM, Xchrom = TRUE,
#'                            Nsim = Nsim, seed = 123)
#' stopifnot(identical(r1.X, r2.X))
#'
#'
#' ### Example 2: Grandparent vs half sib vs uncle
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
#' ### Example 3: X chromosome, granddaughter vs maternal grandfather.
#' y = linearPed(2, sex = c(2, 2))
#' rho = 0.25
#' Nsim = 10
#' estimateTwoLocusIBD(y, c(1,5), rho = rho, Nsim = Nsim, Xchrom = TRUE)
#'
#' # Exact
#' matrix(c(1-rho, rho, rho, 1-rho)/2, ncol = 2)
#'
#' @name estimateIBD
NULL

#' @rdname estimateIBD
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
  map = uniformMap(cM = cM, chromosome = if (Xchrom) 23 else 1)
  
  # Simulate data
  simdata = ibdsim(x, sims = Nsim, ids = ids, map = map, model = "haldane", verbose = verbose, ...)

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

#' @rdname estimateIBD
#' @export
estimateOneLocusIBD = function(x, ids, Nsim, Xchrom = FALSE, verbose = FALSE, ...) {
  st = proc.time()
  
  if (anyNA(match(ids, labels(x))))
    stop2("Unknown ID label: ", setdiff(ids, labels(x)))
  
  # Define map of length 0
  map = uniformMap(cM = 0, chromosome = if (Xchrom) 23 else 1)
  
  # Simulate data
  simdata = ibdsim(x, sims = Nsim, ids = ids, map = map, model = "haldane", verbose = verbose, ...)
  
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


# Utility function: IBD state (0, 1 or 2) for a pair of (non-inbred!) genotypes.
# Each genotype is a pair of alleles.
ibd.state = function(gt1, gt2)
  sum(gt1 %in% gt2)

# Utility function: Jacquard configuration (Sigma 1 - 9) of a pair of genotypes at the same locus.
jacquard.state = function(pat1, mat1, pat2, mat2) {
  if (pat1 == mat1) # Sigma 1,2,3 eller 4
    if (pat2 == mat2) # 1 eller 2
      if (pat1 == pat2) return(1)
      else return(2)
    else
    if (pat1 == pat2 || pat1 == mat2) return(3)
    else return(4)

  if (pat2 == mat2)
    if (pat2 == pat1 || pat2 == mat1) return(5)
    else return(6)

  # If still running: No inbreeding
  ibd = ibd.state(c(pat1, mat1), c(pat2, mat2))
  return((9:7)[ibd + 1])
}


