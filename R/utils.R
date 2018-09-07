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

mergeConsecutiveRows = function(df, mergeBy, segStart=NULL, segEnd=NULL, segLength=NULL) {
  if(nrow(df) < 2) return(df)
  if(length(mergeBy) == 1) 
    mergeBy = df[, mergeBy]
  stopifnot(length(mergeBy) == nrow(df))
      
  runs = rle(mergeBy)
  ends = cumsum(runs$lengths)
  starts = ends - runs$lengths + 1
  
  newdf = df[starts, , drop = FALSE]
  
  if(!is.null(segEnd)) {
    newdf[, segEnd] = df[ends, segEnd]
    if(!is.null(segLength)) {
      newdf[, segLength] = newdf[, segEnd] - newdf[, segStart] 
    }
  }
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
  x[order(x, method="shell")]
}