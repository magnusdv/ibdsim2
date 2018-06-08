# Not exported
meiosis = function(parent, map, model="chi", condition=NULL, skip.recomb=FALSE) { # skip=TRUE returns random strand with no recombination; condition should be NULL or a vector with elements 'locus'(Mb), 'allele' and 'action' (1=force,2=avoid).
  if (condit <- !is.null(condition)) {
    whichStrand = which(condition[["allele"]] == .getAlleles(parent, locus <- condition[["locus"]]))
    startStrand = switch(condition[["action"]],
      switch(length(whichStrand) + 1, stop("Conditional meiosis: Forced allele is not present."), whichStrand, sample.int(2, 1)),
      switch(length(whichStrand) + 1, sample.int(2, 1), 3 - whichStrand, stop("Allele cannot be avoided."))
    )
  } else
    startStrand = sample.int(2, 1)

  if (skip.recomb) return(parent[[startStrand]])
  L.cM = map[nrow(map), "cM"]  # chromosome length in cM

  switch(model,
    haldane = {
      ncross = as.integer(stats::rpois(1, L.cM / 100))
      if (ncross == 0) return(parent[[startStrand]])
      Cx = .sortDouble(stats::runif(ncross, min = 0, max = L.cM))
    },
    chi = {
      m = 4
      nC = stats::rpois(1, L.cM / 50 * (m + 1))    # L.cM/100*2*(m+1); number of potential crossover events
      if (nC == 0) return(parent[[startStrand]])
      C_events = .sortDouble(runif(nC, min = 0, max = L.cM)) # potential crossover positions (N-1 intervals, uniformly distr given nC)
      Cx.bundle = C_events[!as.logical((seq_len(nC) + sample.int(m + 1, 1)) %% (m + 1))]    # Cx events on 4 strand bundle: every (m+1)th
      Cx = Cx.bundle[as.logical(sample.int(2, length(Cx.bundle), replace = T) %% 2)]    # thinning. Each survive with prob=1/2
  })
  cpos = cm2phys(cM_locus = Cx, mapmat = map) # crossover positions
  if (condit)
    startStrand = 2 - (startStrand + sum(cpos < locus)) %% 2  # switches start strand iff sum(cpos < locus) is odd

  recombine(parent[[startStrand]], parent[[3 - startStrand]], cpos)
}


.sortDouble = function(x) x[order(x)] #
