# Not exported
meiosis = function(parent, map, model = "chi", skipRecomb = FALSE) { 
  # skip = TRUE returns random strand with no recombination; 
  
  startStrand = sample.int(2, 1)

  if (skipRecomb) return(parent[[startStrand]])
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
      Cx = Cx.bundle[as.logical(sample.int(2, length(Cx.bundle), replace = TRUE) %% 2)]    # thinning. Each survive with prob=1/2
  })
  
  cpos = cm2phys(cM_locus = Cx, mapmat = map) # crossover positions
  
  # Recombine!
  child = recombine(parent[[startStrand]], parent[[3 - startStrand]], cpos)
  
  # Merge consecutive rows with the same allele
  als = child[, 2]
  if(any(als[-1] == als[-length(als)]))
    child = mergeConsecutiveRows(child, mergeBy = 2)
  
  child
}
