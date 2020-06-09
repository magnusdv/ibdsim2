#' @importFrom pedtools inbreedingLoops descendants 
.CAFs = function(x, sap, loops = NULL) {  # common ancestral founders
  zeros = sap[["0"]]
  ones = c(sap[["1"]], sap[["atleast1"]])
  twos = sap[["2"]]
  caf = fou = setdiff(founders(x, internal = TRUE), zeros)
  fou_one = intersect(fou, ones)
  if (length(fou_one) > 1 || any(twos %in% fou)) 
    return(numeric()) # since founders always have different alleles (by definition)

  fou_desc = lapply(fou, descendants, x = x, internal = TRUE)

  if (length(twos) > 0) {
    .ancestral.founders = function(id) if (id %in% fou) return(id) else fou[sapply(fou_desc, function(ds) id %in% ds)]
    if (is.null(loops)) loops = inbreedingLoops(x)
    if (length(loops) == 0) return(numeric())  # NB: 1) pedigreeLoops uses original ids 2) very inefficient, but fast enough
    bottoms = sapply(loops, "[[", "bottom")
    tops = sapply(loops, "[[", "top")
    for (bot in twos) caf = intersect(caf, unlist(lapply(tops[bottoms == bot], .ancestral.founders)))
  }
  if (length(fou_one) == 1) {
    test_ones = all(ones[ones != fou_one] %in% fou_desc[[match(fou_one, fou)]])
    if (test_ones && fou_one %in% caf) return(fou_one)
    else return(numeric())
  }
  else
    caf[ sapply(caf, function(f) all(ones %in% fou_desc[[match(f, fou)]])) ]  # keep only cafs that have all the 'ones' among their descendants.
}


.pedPaths = function(x, from, to) {
  paths = list()
  lvs = leaves(x, internal = TRUE)

  descend = function(x, from, to, path) {
    if (from == to) {
      paths <- c(paths, list(path))
      return()
    }
    
    if (from %in% leaves) return()
    offs = children(x, from, internal = TRUE)
    for (kid in offs) descend(x, from = kid, to = to, path = c(path, kid))
  }
  descend(x, from, to, path = from)
  paths
}

obligate.carriers = function(x, sap) {
  cafs = .CAFs(x, sap)
  cafPaths = lapply(cafs, function(caf) {
    onePaths = lapply(c(sap[["1"]], sap[["atleast1"]]), function(id) .pedPaths(x, from = caf, to = id))
    twoPaths = lapply(sap[["2"]], function(id) {
      paths = .pedPaths(x, from = caf, to = id)
      temp = list()
      for (i in 1:nrow(loop_pairs <- .comb2(length(paths)))) {
        p1 = paths[[loop_pairs[i, 1]]]; p2 = paths[[loop_pairs[i, 2]]]
        if (p1[length(p1) - 1] != p2[length(p2) - 1]) # true loop only if both parents of 'two' are included
          temp = c(temp, list(c(p1, p2)))
      }
      temp
    }) # twoPaths for this CAF is a list of |two| elements: For each two-ID the list of possible loops from CAF to ID
    allcomb = expand.grid(c(twoPaths, onePaths)) # each row here consists of one of each
    lapply(1:nrow(allcomb), function(i) sort.default(unique(unlist(allcomb[i, ]))))
  })
  obligs = unique(unlist(cafPaths, recursive = FALSE))
  obligs = lapply(obligs, setdiff, c(sap[["atleast1"]], sap[["2"]]))
  obligs = obligs[sapply(obligs, function(vec) !any(sap[["0"]] %in% vec))] # and paths containing 0-indivs
  if (length(obligs) == 0) 
    stop2("Cannot find any possible sets of obligate carriers")
  
  keep = lapply(seq_len(length(obligs)), function(i) {
    for (other in obligs[-i]) 
      if (all(other %in% obligs[[i]])) return(FALSE)
    return(TRUE)
  })
  obligs[unlist(keep)]
}
