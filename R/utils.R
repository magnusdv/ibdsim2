stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}

# Test that input is a positive (or similar) integer.
is_count = function(x, minimum = 1) {
  isTRUE(length(x) == 1 &&
           (is.integer(x) || (is.numeric(x) && x == as.integer(x))) &&
           x >= minimum)
}

mergeConsecutiveRows = function(df, mergeBy) {
  if(nrow(df) < 2) return(df)
  if(length(mergeBy) == 1) 
    mergeBy = df[, mergeBy]
  stopifnot(length(mergeBy) == nrow(df))
      
  runs = rle(mergeBy)
  ends = cumsum(runs$lengths)
  starts = ends - runs$lengths + 1
  
  df[starts, , drop = FALSE]
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

.getAlleles = function(chromdata, posvec) {
  posvec[posvec < 0] = 0
  rbind(pos2allele(chromdata[[1]], posvec),
    pos2allele(chromdata[[2]], posvec))
}


pos2allele = function(haplo, posvec) { # haplo = matrix with 2 columns (breaks - allele)
  indices = findInterval(posvec, haplo[, 1])
  haplo[indices, 2]
}

# Unused at the moment
#locus2marker = function(x, h, locus) {
#  marker = t.default(sapply(h, .getAlleles, posvec = locus))
#  paramlink::setMarkers(x, marker)
#}

.printSAP = function(sap) {
  if (!is.null(two <- sap[["2"]])) cat("  Two copies:", paste(two, collapse = ", "), "\n")
  if (!is.null(atl1 <- sap[["atleast1"]])) cat("  At least one copy:", paste(atl1, collapse = ", "), "\n")
  if (!is.null(one <- sap[["1"]])) cat("  One copy:", paste(one, collapse = ", "), "\n")
  if (!is.null(atm1 <- sap[["atmost1"]])) cat("  At most one copy:", paste(atm1, collapse = ", "), "\n")
  if (!is.null(zero <- sap[["0"]])) cat("  Zero copies:", paste(zero, collapse = ", "), "\n")
  cat("\n")
}

.sortDouble = function(x) {
  x[order(x, method = "shell")]
}

.comb2 = function(n) {
  if (n < 2) return(matrix(nrow = 0, ncol = 2))
  v1 = rep.int(seq_len(n - 1), (n - 1):1)
  v2 = NULL
  for (i in 2:n) v2 = c(v2, i:n)
  cbind(v1, v2, deparse.level = 0)
}

`%||%` = function(x, y) {
  if (is.null(x)) y
  else x
}

vecToRanges = function(x) {
  xInt = suppressWarnings(as.integer(x))
  if(any(is.na(xInt) | x != xInt))
    return(x)

  x = xInt
  conseqs = unname(split(x, cumsum(c(0, diff(x) != 1))))
  rngs = sapply(conseqs, function(v) {
    if(length(v) == 1) as.character(v)
    else if(length(v) == 2) toString(v)
    else sprintf("%d-%d", min(v), max(v))
  })
  paste(rngs, collapse = ", ")
}