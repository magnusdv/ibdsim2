
alleleFlow = function(x, ids, addState = TRUE) {
  if(!inherits(x, "genomeSim"))
    stop2("Argument `x` must be a `genomeSim` object. Received: ", class(x))
  
  xids = extractIds(x)
  if(!all(ids %in% xids))
    stop2("Unknown ID label: ", setdiff(ids, xids))
  
  n = length(ids)
  
  # If nothing to do, return early
  if(setequal(ids, xids)) {
    if(n > 2 || !addState) 
      return(x)
    if(n == 1 && "Aut" %in% colnames(x))
      return(x)
    if(n == 2 && "Sigma" %in% colnames(x))
      return(x)
  }
  
  # Allele columns of selected ids
  cols = paste(rep(ids, each = 2), c("p", "m"), sep = ":")
  y = mergeSegments(x, by = cols)
  
  # Keep only allele columns of selected ids
  y = y[, c(1:4, which(colnames(x) %in% cols)), drop = FALSE]
  
  if(addState)
    y = addStates(y)
  
  y
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
      else
        byEq = rowSums(x[-1, by, drop = FALSE] == x[-k, by, drop = FALSE]) == byLen
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
    adjac = x[-1, 'start'] == x[-k, 'end']
    mergeRow = mergeRow & adjac
  }
 
  if(!any(mergeRow))
    return(x)
  
  fromRow = which(c(TRUE, !mergeRow))
  toRow = c(fromRow[-1] - 1, k)
  
  y = x[fromRow, , drop = FALSE]
  y[, 'end'] = x[toRow, 'end']
  y[, 'length'] = y[, 'end'] - y[, 'start'] 
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
#' segmentStats(segs)
#'
#' @importFrom stats quantile
#' @export
segmentStats = function(x, quantiles = c(0.025, 0.5, 0.975), returnAll = FALSE) {
  
  if(is.matrix(x))
    x = list(x)
  
  # List of lengths
  lenDat = lapply(x, function(m) m[, 'length'])
  
  # Collect data for each sim
  perSim = data.frame(
    `Count`    = lengths(lenDat),
    `Total`    = vapply(lenDat, sum, FUN.VALUE = numeric(1)), 
    `Average`  = vapply(lenDat, safeMean, FUN.VALUE = numeric(1)), 
    `Shortest` = vapply(lenDat, safeMin, FUN.VALUE = numeric(1)),
    `Longest`  = vapply(lenDat, safeMax, FUN.VALUE = numeric(1)))
  
  # Summarising function
  sumfun = function(v) c(mean = safeMean(v), sd = sd(v), min = safeMin(v), quantile(v, quantiles), max = safeMax(v))
  
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


# TODO: Remove. Superseded by the new `mergeSegments()`
mergeConsecutiveSegments = function(df, mergeBy, segStart = "start", 
                                    segEnd = "end", segLength = "length") {
  N = nrow(df)
  if(N < 2) return(df)
  
  if(length(mergeBy) == 1)
    mergeBy = df[, mergeBy]
  else if(length(mergeBy) != N)
    mergeBy = apply(df[, mergeBy], 1, paste, collapse = ":")
  
  newSeg = c(TRUE, df[-N, segEnd] != df[-1, segStart])
  mergeBy = paste(mergeBy, cumsum(newSeg), sep = "_chunk")
  
  # Runs of consecutive segs; row numbers of start/end
  runs = rle(mergeBy)
  ends = cumsum(runs$lengths)
  starts = ends - runs$lengths + 1
  
  # Keep only rows starting a new segment
  newdf = df[starts, , drop = FALSE]
  
  # Fix endpoints
  newdf[, segEnd] = df[ends, segEnd]
  
  # Fix lengths if present
  if(!is.null(segLength) && segLength %in% colnames(df))
    newdf[, segLength] = newdf[, segEnd] - newdf[, segStart] 
  
  newdf
}



# Add columns with IBD states (if 1 or 2 ids)
addStates = function(x) {
  
  nids = (ncol(x) - 4)/2 
  Xchrom = any(x[, 1] == 23)
  
  if(nids == 1) {
    a1 = x[, 5]; a2 = x[, 6]
    aut = a1 == a2 
    if(Xchrom) # convention: hemizygous = autozygous
      aut = aut | a1 == 0 
    aut = as.numeric(aut)  # x is a numeric (not integer) matrix
    x = cbind(x, Aut = aut)
  }
  else if(nids == 2) {
    a1 = x[, 5]; a2 = x[, 6] 
    b1 = x[, 7]; b2 = x[, 8]
    
    IBD = ibdState(a1, a2, b1, b2, Xchrom = Xchrom)
    Sigma = jacquardState(a1, a2, b1, b2, Xchrom = Xchrom)
    
    x = cbind(x, IBD = IBD, Sigma = Sigma)
  }
  
  x
}

# Return the condensed identity ("jacquard") state given 4 alleles:
# a1 and a2 from individual 1; b1 and b2 from individual 2
jacquardState = function(a1, a2, b1, b2, Xchrom = FALSE) { 
  
  # Identity states on X are follow the autosomal ordering,
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

