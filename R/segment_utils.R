
alleleFlow = function(x, ids, addState = TRUE) {
  if(!inherits(x, "genomeSim"))
    stop2("Argument `x` must be a `genomeSim` object. Received: ", class(x))
  
  xids = extractIds(x)
  if(!all(ids %in% xids))
    stop2("Unknown ID label: ", setdiff(ids, xids))
  
  n = length(ids)
  clnms = colnames(x)
  
  # If nothing to do, return early
  if(setequal(ids, xids)) {
    if(n > 2 || !addState) 
      return(x)
    if(n == 1 && "Aut" %in% clnms)
      return(x)
    if(n == 2 && "Sigma" %in% clnms)
      return(x)
  }
  
  # Allele columns of selected ids
  cols = paste(rep(ids, each = 2), c("p", "m"), sep = ":")
  y = mergeSegments(x, by = cols)
  
  # Keep only allele columns of selected ids
  y = y[, c(clnms[1:5], cols), drop = FALSE]
  
  if(addState)
    y = addStates(y, cols)
  
  y
}


stackIntervals = function(start, end, chrom = NULL, sort = TRUE) {
  
  n = length(start)
  if(n < 2)
    return(rep_len(1, n))
  
  lay = integer(n)
  
  if(!is.null(chrom)) {
    for(chr in unique.default(chrom)) {
      idx = chrom == chr
      lay[idx] = stackIntervals(start = start[idx], end = end[idx], sort = sort)
    }
    return(lay)
  }
    
  if(any(start > end))
    stop2("Inverted segment")
  
  if(sort) {
    ord = order(start, -end, method = "shell")
    start = start[ord]
    end = end[ord]
  }
  
  # Current end of each layer
  layerEnds = numeric(n) # theoretical max
  layerEnds[1] = end[1]
  
  lay[1] = 1
  for (k in seq.int(2, n)) {
    for(i in 1:n) {
      if(layerEnds[i] < start[k]) {
        lay[k] = i
        layerEnds[i] = end[k]
        break
      }
    }
  }
  
  if(sort)
    lay = lay[order(ord, method = "shell")]
  lay
}

groupOverlaps = function(x) { # assumes sorted on chrom, start, end
  
  # Storage vector for group number and group size
  g = n = gidx = numeric(nrow(x))
  chrom = x[, 1]
  start = x[, 2]
  end = x[, 3]
  
  for(chr in unique.default(chrom)) {
    idx = chrom == chr
    sta = start[idx]
    sto = end[idx]
    gchr = cumsum(cummax(c(0, sto[-length(sto)])) <= sta)
    tb = tabulate(gchr)
    g[idx] = gchr
    n[idx] = tb[gchr]
    gidx[idx] = unlist(lapply(tb, seq_len), use.names = FALSE)
  }
  
  cbind(x, group = paste(chrom, g, sep = "-"), gsize = n, gidx = gidx)
}

# Merge segments (on the same chrom) with equal entries in `by` (single vector or column names)
mergeSegments = function(x, by = NULL, checkAdjacency = FALSE) {
  k = nrow(x)
  if(k < 2)
    return(x)
  
  # Rows where chromosome is unchanged
  chrEq = x[-1, 'chrom'] == x[-k, 'chrom']
  mergeRow = chrEq
  
  # Rows where `by` elements don't change
  if(!is.null(by)) {
    byLen = length(by)
    
    if(is.character(by) && all(by %in% colnames(x))) { # Interpret as column names
      if(byLen == 1)
        byEq = x[-1, by] == x[-k, by]
      else {
        chck = x[-1, by, drop = FALSE] == x[-k, by, drop = FALSE]
        byEq = .rowSums(chck, m = k - 1, n = byLen) == byLen
      }
    }
    else if(byLen == k && (is.numeric(by) || is.character(by))) {
      byEq = by[-1] == by[-k]
    }
    else
      stop2("Wrong format or length of argument `by`")
    
    mergeRow = chrEq & !is.na(byEq) & byEq
  }
 
  # Segment adjacency
  if(checkAdjacency) {
    adjac = x[-1, 'startMB'] == x[-k, 'endMB']
    mergeRow = mergeRow & adjac
  }
 
  if(!any(mergeRow))
    return(x)
  
  fromRow = which(c(TRUE, !mergeRow))
  toRow = c(fromRow[-1] - 1, k)
  
  y = x[fromRow, , drop = FALSE]
  y[, 'endMB'] = x[toRow, 'endMB']
  y[, 'endCM'] = x[toRow, 'endCM']
  
  rownames(y) = NULL
  y
}


#' Summary statistics for identified segments
#'
#' Compute summary statistics for segments identified by [findPattern()].
#'
#' @param x A list of matrices produced with [findPattern()].
#' @param quantiles A vector of quantiles to include in the summary.
#' @param returnAll A logical, by default FALSE. If TRUE, the output includes a
#'   vector `allSegs` containing the lengths of all segments in all simulations.
#' @param unit Either "mb" (megabases) or "cm" (centiMorgan); the length unit
#'   for genomic segments.
#'   
#' @return A list containing a data frame `perSim`, a matrix `summary` and (if
#'   `returnAll` is TRUE) a vector `allSegs`.
#'
#'   Variables used in the output:
#'
#'   * `Count`: The total number of segments in a simulation
#'
#'   * `Total`: The total sum of the segment lengths in a simulation
#'
#'   * `Average`: The average segment lengths in a simulation
#'
#'   * `Shortest`: The length of the shortest segment in a simulation
#'
#'   * `Longest`: The length of the longest segment in a simulation
#'
#'   * `Overall` (only in `summary`): A summary of all segments from all
#'   simulations
#'
#' @seealso [findPattern()]
#'
#' @examples
#' x = nuclearPed(3)
#' sims = ibdsim(x, N = 2, map = uniformMap(M = 2), model = "haldane", seed = 1729)
#'
#' # Segments where all siblings carry the same allele
#' segs = findPattern(sims, pattern = list(carriers = 3:5))
#'
#' # Summarise
#' segmentStats(segs, unit = "mb")
#' 
#' # The unit does not matter in this case (since the map is trivial)
#' segmentStats(segs, unit = "cm")
#'
#' @importFrom stats quantile
#' @export
segmentStats = function(x, quantiles = c(0.025, 0.5, 0.975), returnAll = FALSE, unit = "mb") {
  
  if(is.matrix(x))
    x = list(x)
  
  # Names of start/end columns
  startCol = switch(unit, mb = "startMB", cm = "startCM")
  endCol = switch(unit, mb = "endMB", cm = "endCM")
  
  # List of lengths
  lenDat = lapply(x, function(s) s[, endCol] - s[, startCol])
  
  # Collect data for each sim
  perSim = data.frame(
    `Count`    = lengths(lenDat),
    `Total`    = vapply(lenDat, sum, FUN.VALUE = numeric(1)), 
    `Average`  = vapply(lenDat, safeMean, FUN.VALUE = numeric(1)), 
    `Shortest` = vapply(lenDat, safeMin, FUN.VALUE = numeric(1)),
    `Longest`  = vapply(lenDat, safeMax, FUN.VALUE = numeric(1)))
  
  # Summarising function
  sumfun = function(v) c(mean = safeMean(v), sd = sd(v), min = safeMin(v), 
                         quantile(v, quantiles), max = safeMax(v))
  # Summary
  sumList = lapply(perSim, sumfun)
  
  # Add stats of overall lengths
  allSegs = unlist(lenDat, use.names = FALSE)
  sumList$Overall = sumfun(allSegs)
  
  # Return as data frames
  res = list(perSim = perSim, summary = do.call(rbind, sumList))
  if(returnAll)
    res$allSegs = allSegs
  
  res
}


# Add columns with IBD states (if 1 or 2 ids)
addStates = function(x, acols) {
  
  nc = length(acols)
  Xchrom = any(x[, "chrom"] == 23)
  
  if(nc == 2) {
    a1 = x[, acols[1]]
    a2 = x[, acols[2]]
    aut = a1 == a2 
    if(Xchrom) # convention: hemizygous = autozygous
      aut = aut | a1 == 0 
    aut = as.numeric(aut)  # x is a numeric (not integer) matrix
    x = cbind(x, Aut = aut)
  }
  else if(nc == 4) {
    a1 = x[, acols[1]]
    a2 = x[, acols[2]]
    b1 = x[, acols[3]]
    b2 = x[, acols[4]]
    
    IBD = ibdState(a1, a2, b1, b2, Xchrom = Xchrom)
    Sigma = jacquardState(a1, a2, b1, b2, Xchrom = Xchrom)
    
    x = cbind(x, IBD = IBD, Sigma = Sigma)
  }
  else
    stop2("addStates requires 2 or 4 columns, not: ", nc)
  
  x
}

# Return the condensed identity ("jacquard") state given 4 alleles:
# a1 and a2 from individual 1; b1 and b2 from individual 2
jacquardState = function(a1, a2, b1, b2, Xchrom = FALSE) { 
  
  # Identity states on X follow the autosomal ordering,
  # after replacing hemizygous alleles with autozygous ones. 
  # See https://github.com/magnusdv/ribd
  if(Xchrom) {
    azero = a1 == 0
    bzero = b1 == 0
    a1[azero] = a2[azero]
    b1[bzero] = b2[bzero]
  }
  
  inb1 = a1 == a2
  inb2 = b1 == b2
  pa = a1 == b1
  ma = a2 == b2
  d1 = a1 == b2
  d2 = a2 == b1
  
  horiz = inb1 + inb2
  verts = pa + ma + d1 + d2
  
  res = integer(length(a1))
  
  res[horiz == 2 & pa] = 1
  res[horiz == 2 & !pa] = 2
  
  res[inb1 & !inb2 & (pa | ma)] = 3
  res[inb1 & !inb2 & verts == 0] = 4
  
  res[!inb1 & inb2 & (pa | ma)] = 5
  res[!inb1 & inb2 & verts == 0] = 6
  
  res[horiz == 0 & verts == 2] = 7
  res[verts == 1] = 8
  res[horiz + verts == 0] = 9
  res
}

ibdState = function(a1, a2, b1, b2, Xchrom = FALSE) { 
  
  inb1 = a1 == a2
  inb2 = b1 == b2
  pa = a1 == b1
  ma = a2 == b2
  d1 = a1 == b2
  d2 = a2 == b1
  
  if(Xchrom) {
    hem1 = a1 == 0
    hem2 = b1 == 0
  
    pa = pa & !hem1 & !hem2
    d1 = d1 & !hem1
    d2 = d2 & !hem2
  }
  
  res = as.integer(pa + ma + d1 + d2)
  res[inb1 | inb2] = NA_integer_
  res
}

