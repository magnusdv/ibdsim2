twoLocusIBD = function(x, ind1, ind2, rho=NULL, cM=NULL, Nsim, Xchrom=FALSE, verbose=TRUE, ...) {
  if ((err <- ind1) > x$nInd || (err <- ind2) > x$nInd) stop(paste("The pedigree has no individual with label", err))
  if (is.null(cM) + is.null(rho) != 1) stop("Exactly one of the parameters 'cM' and 'rho' must be NULL.")
  starttime = proc.time()

  if (is.null(cM)) {
    if (rho < 0 | rho > 0.5) stop("Recombination rate 'rho' is outside of interval [0, 0.5].")
    cM = -50 * log(1 - 2 * rho)
  }
  if (verbose) cat("Locus distance:", cM, "centiMorgan\n")

  if (cM == Inf) {
    if (verbose) cat("Analysis of unlinked loci.\n")
    m1 = oneLocusIBD(x, ind1, ind2, Nsim = Nsim, Xchrom = Xchrom, verbose = verbose, ...)
    m2 = oneLocusIBD(x, ind1, ind2, Nsim = Nsim, Xchrom = Xchrom, ...) # dont repeat verbose output
    res = outer(m1, m2)
    return(res)
  }

  map = uniformMap(cM = cM, chromosome = if (Xchrom) 23 else 1)
  simdata = ibdsim(x, map = map, sims = Nsim, model = "haldane", verbose = verbose, ...)

  # Setup for X chromosomal analysis
  if (Xchrom) {
    sex1 = x$pedigree[ind1, "SEX"] # NB: ibdsim raises error if not labels are 1,2,3,...
    sex2 = x$pedigree[ind2, "SEX"]

    if (sex1 == 2 && sex2 == 1) { # reverse indiv order in this case (simplifies code below)
      inds = c(ind1, ind2)
      ind1 = inds[2]; ind2 = inds[1]
      sex1 = 1; sex2 = 2
    }
  }

  # Define appropriate function for analysing a single simulation: Output is a pair (IBD_locus1, IBD_locus2)
  if (Xchrom && sex1 == 1 && sex2 == 1)
    single.sim.ibd = function(sim) {
      i1.mat = sim[[c(1, ind1, 2)]] # maternal chromosome of indiv 1 (matrix with 2 columns: breakpoint position; allele)
      i2.mat = sim[[c(1, ind2, 2)]]

      # locus 1
      m1_ibd = i1.mat[1, 2] == i2.mat[1, 2] # first row (first locus), second column (allele)

      # locus 2
      m2_ibd =  i1.mat[dim(i1.mat)[1], 2] == i2.mat[dim(i2.mat)[1], 2] # last row (second locus), second column (allele)
      c(m1_ibd, m2_ibd)
    }
  else if (Xchrom && sex1 == 1) # includes sex1=2, sex2=1, because of swapping above
    single.sim.ibd = function(sim) {
      i1.mat = sim[[c(1, ind1, 2)]] # maternal chromosome of indiv 1 (matrix with 2 columns: breakpoint position; allele)
      i2.pat = sim[[c(1, ind2, 1)]]
      i2.mat = sim[[c(1, ind2, 2)]]

      # locus 1
      ind1.mat.m1 = i1.mat[1, 2]
      ind2.pat.m1 = i2.pat[1, 2]
      ind2.mat.m1 = i2.mat[1, 2]
      m1_ibd =  ind1.mat.m1 %in% c(ind2.pat.m1, ind2.mat.m1)

      # locus 2
      ind1.mat.m2 = i1.mat[dim(i1.mat)[1], 2] # last row (second locus), second column (allele)
      ind2.pat.m2 = i2.pat[dim(i2.pat)[1], 2]
      ind2.mat.m2 = i2.mat[dim(i2.mat)[1], 2]
      m2_ibd =  ind1.mat.m2 %in% c(ind2.pat.m2, ind2.mat.m2)

      c(m1_ibd, m2_ibd)
    }
  else
    single.sim.ibd = function(sim) {
      i1.pat = sim[[c(1, ind1, 1)]] # paternal chromosome of indiv 1 (matrix with 2 columns: breakpoint position; allele)
      i1.mat = sim[[c(1, ind1, 2)]]
      i2.pat = sim[[c(1, ind2, 1)]]
      i2.mat = sim[[c(1, ind2, 2)]]

      # locus 1
      ind1.pat.m1 = i1.pat[1, 2]
      ind1.mat.m1 = i1.mat[1, 2]
      ind2.pat.m1 = i2.pat[1, 2]
      ind2.mat.m1 = i2.mat[1, 2]
      m1_ibd =  ibd.status(c(ind1.pat.m1, ind1.mat.m1), c(ind2.pat.m1, ind2.mat.m1))

      # locus 2
      ind1.pat.m2 = i1.pat[dim(i1.pat)[1], 2]
      ind1.mat.m2 = i1.mat[dim(i1.mat)[1], 2]
      ind2.pat.m2 = i2.pat[dim(i2.pat)[1], 2]
      ind2.mat.m2 = i2.mat[dim(i2.mat)[1], 2]
      m2_ibd =  ibd.status(c(ind1.pat.m2, ind1.mat.m2), c(ind2.pat.m2, ind2.mat.m2))

      c(m1_ibd, m2_ibd)
    }


  # Analyse simulations
  allsims.ibd = vapply(simdata, single.sim.ibd, numeric(2))

  # Frequency table of results
  max.ibd = if (Xchrom && sex1 == 1) 1 else 2
  res = matrix(ncol = max.ibd + 1, nrow = max.ibd + 1, dimnames = rep(list(paste0("ibd", 0:max.ibd)), 2))
  for (i in 0:max.ibd) for (j in 0:max.ibd)
      res[i + 1, j + 1] = sum((allsims.ibd[1, ] == i) & (allsims.ibd[2, ] == j))

  if (verbose) cat("IBD analysis finished. Total time used:", (proc.time() - starttime)[["elapsed"]], "seconds.\n")
  res / Nsim
}


oneLocusIBD = function(x, ind1, ind2, Nsim, Xchrom=FALSE, verbose=TRUE, ...) {
  if ((err <- ind1) > x$nInd || (err <- ind2) > x$nInd) stop(paste("The pedigree has no individual with label", err))

  # Setup for X chromosomal analysis
  if (Xchrom) {
    sex1 = x$pedigree[ind1, "SEX"] # NB: ibdsim raises error if not labels are 1,2,3,...
    sex2 = x$pedigree[ind2, "SEX"]

    if (sex1 == 2 && sex2 == 1) { # reverse indiv order in this case (simplifies code below)
      inds = c(ind1, ind2)
      ind1 = inds[2]; ind2 = inds[1]
      sex1 = 1; sex2 = 2
    }
  }

  # Define appropriate function for analysing a single simulation: Output is a single integer (0,1,2)
  if (Xchrom && sex1 == 1 && sex2 == 1)
    single.sim.ibd = function(sim) {
      i1.mat = sim[[c(1, ind1, 2)]] # maternal chromosome of indiv 1. 1*2 matrix: [0, allele]
      i2.mat = sim[[c(1, ind2, 2)]]
      i1.mat[1, 2] == i2.mat[1, 2] # TRUE iff IBD = 1
    }
  else if (Xchrom && sex1 == 1) # includes sex1=2, sex2=1, because of swapping above
    single.sim.ibd = function(sim) {
      i1.mat = sim[[c(1, ind1, 2)]] # maternal chromosome of indiv 1
      i2.pat = sim[[c(1, ind2, 1)]]
      i2.mat = sim[[c(1, ind2, 2)]]
      # each of the above is a 1*2 matrix: [0, allele]

      i1.mat[1, 2] %in% c(i2.pat[1, 2], i2.mat[1, 2])
    }
  else
    single.sim.ibd = function(sim) {
      i1.pat = sim[[c(1, ind1, 1)]] # paternal chromosome of indiv 1 (matrix with 2 columns: breakpoint position; allele)
      i1.mat = sim[[c(1, ind1, 2)]]
      i2.pat = sim[[c(1, ind2, 1)]]
      i2.mat = sim[[c(1, ind2, 2)]]
      # each of the above is a 1*2 matrix: [0, allele]

      ibd.status(c(i1.pat[1, 2], i1.mat[1, 2]), c(i2.pat[1, 2], i2.mat[1, 2]))
    }

  starttime = proc.time()

  map = uniformMap(cM = 0, chromosome = if (Xchrom) 23 else 1)
  simdata = ibdsim(x, map = map, sims = Nsim, model = "haldane", verbose = verbose, ...)

  # Analyse simulations
  allsims.ibd = vapply(simdata, single.sim.ibd, numeric(1))

  ibdlevels = if (Xchrom && sex1 == 1) 0:1 else 0:2
  res = table(factor(allsims.ibd, levels = ibdlevels, labels = paste0("ibd", ibdlevels)))
  res = c(res) # strips table attributes, but keeps names. (Hadley: Bad style.)
  if (verbose) cat("IBD analysis finished. Total time used:", (proc.time() - starttime)[["elapsed"]], "seconds.\n")
  res / Nsim
}


# Utility function: IBD status (0, 1 or 2) for a pair of (non-inbred!) genotypes.
# Each genotype is a pair of alleles.
ibd.status = function(gt1, gt2)
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
  ibd = ibd.status(c(pat1, mat1), c(pat2, mat2))
  return((9:7)[ibd + 1])
}

# Observed frequencies of two-locus Jacquard condensed identity states.
# Similar to twoLocusIBD bortsett, except that jacquard.state() is used in place of ibd.status(), and the use of table() at the end.

twoLocusJacquard = function(x, ind1, ind2, rho=NULL, cM=NULL, Nsim, verbose=TRUE, ...) {
  if (is.null(cM) + is.null(rho) != 1) stop("Exactly one of the parameters 'cM' and 'rho' must be NULL.")
  starttime = proc.time()

  if (is.null(cM)) {
    if (rho < 0 | rho > 0.5) stop("Recombination rate 'rho' is outside of interval [0, 0.5].")
    cM = -50 * log(1 - 2 * rho)
  }
  if (verbose) cat("Locus distance:", cM, "centiMorgan\n")

  if (cM == Inf) {
    if (verbose) cat("Analysing of unlinked loci.\n")
    m1 = oneLocusJacquard(x, ind1, ind2, Nsim = Nsim, verbose = verbose, ...)
    m2 = oneLocusJacquard(x, ind1, ind2, Nsim = Nsim, ...) # dont repeat verbose output
    res = outer(m1, m2)
    return(res)
  }

  map = uniformMap(cM = cM)
  simdata = ibdsim(x, sims = Nsim, map = map, model = "haldane", ...)

  sim.jacquard = sapply(simdata, function(sim) {
    i1.pat = sim[[c(1, ind1, 1)]] # paternal chromosome of indiv 1 (matrix with 2 columns: breakpoint position; allele)
    i1.mat = sim[[c(1, ind1, 2)]]
    i2.pat = sim[[c(1, ind2, 1)]]
    i2.mat = sim[[c(1, ind2, 2)]]

    # locus 1
    ind1.pat.m1 = i1.pat[1, 2]
    ind1.mat.m1 = i1.mat[1, 2]
    ind2.pat.m1 = i2.pat[1, 2]
    ind2.mat.m1 = i2.mat[1, 2]
    m1_Sigma =  jacquard.state(ind1.pat.m1, ind1.mat.m1, ind2.pat.m1, ind2.mat.m1)

    # locus 2
    ind1.pat.m2 = i1.pat[dim(i1.pat)[1], 2]
    ind1.mat.m2 = i1.mat[dim(i1.mat)[1], 2]
    ind2.pat.m2 = i2.pat[dim(i2.pat)[1], 2]
    ind2.mat.m2 = i2.mat[dim(i2.mat)[1], 2]
    m2_Sigma =  jacquard.state(ind1.pat.m2, ind1.mat.m2, ind2.pat.m2, ind2.mat.m2)

    c(m1_Sigma, m2_Sigma)
  })

  res = table(factor(sim.jacquard[1, ], levels = 1:9), factor(sim.jacquard[2, ], levels = 1:9))
  if (verbose) cat("Jacquard analysis finished.\nTotal time used:", (proc.time() - starttime)[["elapsed"]], "seconds.\n---------------\n")
  as.matrix(res) / Nsim
}


oneLocusJacquard = function(x, ind1, ind2, Nsim, verbose=TRUE, ...) {
  if ((err <- ind1) > x$nInd || (err <- ind2) > x$nInd) stop(paste("The pedigree has no individual with label", err))

  sim.jacquard = function(sim) {
    i1.pat = sim[[c(1, ind1, 1)]] # paternal chromosome of indiv 1 (a 1*2 matrix: [0, allele])
    i1.mat = sim[[c(1, ind1, 2)]]
    i2.pat = sim[[c(1, ind2, 1)]]
    i2.mat = sim[[c(1, ind2, 2)]]

    jacquard.state(i1.pat[1, 2], i1.mat[1, 2], i2.pat[1, 2], i2.mat[1, 2])
  }

  starttime = proc.time()

  map = uniformMap(cM = 0)
  simdata = ibdsim(x, map = map, sims = Nsim, model = "haldane", verbose = verbose, ...)

  # Analyse simulations
  allsims.jacq = vapply(simdata, sim.jacquard, numeric(1))

  # Output
  res = table(factor(allsims.jacq, levels = 1:9))
  res = c(res) # strips table attributes, but keeps names. (Hadley: Bad style.)
  if (verbose) cat("Jacquard analysis finished. Total time used:", (proc.time() - starttime)[["elapsed"]], "seconds.\n")
  res / Nsim
}
