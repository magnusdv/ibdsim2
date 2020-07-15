#' @importFrom stats rpois runif
meiosis = function(parent, map, model = "chi", skipRecomb = FALSE) { 
  startStrand = sample.int(2, 1)

  if (skipRecomb) # return random strand with no recombination
    return(parent[[startStrand]])
  
  L.cM = map[nrow(map), "cM"]  # chromosome length in cM

  switch(model,
    haldane = {
      ncross = rpois(1, L.cM / 100)
      if (ncross == 0) 
        return(parent[[startStrand]])
      Cx = .sortDouble(runif(ncross, min = 0, max = L.cM))
    },
    chi = {
      m = 4
      nC = rpois(1, L.cM / 50 * (m + 1))    # L.cM/100*2*(m+1); number of potential crossover events
      if (nC == 0) 
        return(parent[[startStrand]])
      C_events = .sortDouble(runif(nC, min = 0, max = L.cM)) # potential crossover positions (N-1 intervals, uniformly distr given nC)
      Cx.bundle = C_events[!as.logical((seq_len(nC) + sample.int(m + 1, 1)) %% (m + 1))]    # Cx events on 4 strand bundle: every (m+1)th
      Cx = Cx.bundle[as.logical(sample.int(2, length(Cx.bundle), replace = TRUE) %% 2)]    # thinning. Each survive with prob=1/2
  })
  
  cpos = convertPos(cM = Cx, map = map) # crossover positions
  
  # Recombine!
  child = recombine(parent[[startStrand]], parent[[3 - startStrand]], cpos)
  
  # Merge consecutive rows with the same allele (needed in inbred peds)
  als = child[, 2]
  nr = length(als)
  if(nr > 1 && any(als[-1] == als[-nr])) {
    runs = rle(als)
    ends = cumsum(runs$lengths)
    starts = ends - runs$lengths + 1

    child = child[starts, , drop = FALSE]
  }
  
  child
}

