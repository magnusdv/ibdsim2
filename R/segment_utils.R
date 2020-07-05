segmentSummary = function(x, ids, addState = TRUE) {
  if(!inherits(x, "genomeSim"))
    stop2("Argument `x` must be a `genomeSim` object. Received: ", class(x))
  
  xids = extractIdsFromSegmentSummary(x)
  if(!all(ids %in% xids))
    stop2("Unknown ID label: ", setdiff(ids, xids))
  
  n = length(ids)
  
  if(setequal(ids, xids)) {
    if(n > 2 || !addState) 
      return(x)
    if(n == 1 && "Aut" %in% colnames(x))
      return(x)
    if(n == 2 && "Sigma" %in% colnames(x))
      return(x)
  }
  
  # Allele columns of selected ids
  colnms = paste(rep(ids, each = 2), c("p", "m"), sep = ":")
  
  # Merge identical rows
  y = mergeAdjacent(x, vec = apply(x[, colnms], 1, paste, collapse = " "))
  y = y[, c(1:4, match(colnms, colnames(x))), drop = FALSE]
  
  if(addState)
    y = addStates(y)
  
  y
}

# Merge adjacent segments with equal `vec` entry and equal chrom
mergeAdjacent = function(x, vec) {
  k = nrow(x)
  
  if(length(vec) == 1 && is.character(vec))
    vec = x[, vec]
  else if(length(vec) != k)
    stop2("Incompatible input")
  
  chr = x[, 'chrom']
  
  vecEq = vec[-1] == vec[-k]
  chrEq = chr[-1] == chr[-k]
  
  mergeRow = vecEq & chrEq
  if(!any(mergeRow))
    return(x)
  
  fromRow = which(c(TRUE, !mergeRow))
  toRow = c(fromRow[-1] - 1, k)
  
  
  y = x[fromRow, , drop = FALSE]
  y[, 'end'] = x[toRow, 'end']
  y[, 'length'] = y[, 'end'] - y[, 'start'] 
  y
}


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
  if(!is.null(segLength) && segLength %in% names(df))
    newdf[, segLength] = newdf[, segEnd] - newdf[, segStart] 
  
  newdf
}



# Add columns with IBD states (if 1 or 2 ids)
addStates = function(x) {
  if(ncol(x) == 8) { # Pairwise
    # Jacquard state
    Sigma = jacquardState(a1 = x[, 5], a2 = x[, 6], b1 = x[, 7], b2 = x[, 8])
    
    # IBD state (0 <-> 9; 1 <-> 8; 2 <-> 7; otherwise NA)
    IBD = 9 - Sigma
    IBD[IBD > 2] = NA
    
    return(cbind(x, IBD = IBD, Sigma = Sigma))
  } 
  else if(ncol(x) == 6) {
    # Autozygosity state
    a = as.numeric(x[, 5] == x[, 6]) # since x is numeric
    return(cbind(x, Aut = a))
  }
  else 
    return(x)
}

# Return the condensed identity ("jacquard") state given 4 alleles:
# a1 and a2 from individual 1; b1 and b2 from individual 2
jacquardState = function(a1, a2, b1, b2) { 
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


