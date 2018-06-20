
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
locus2marker = function(x, h, locus) {
  marker = t.default(sapply(h, .getAlleles, posvec = locus))
  paramlink::setMarkers(x, marker)
}

.printSAP = function(sap) {
  if (!is.null(two <- sap[["2"]])) cat("  Two copies:", paste(two, collapse = ", "), "\n")
  if (!is.null(atl1 <- sap[["atleast1"]])) cat("  At least one copy:", paste(atl1, collapse = ", "), "\n")
  if (!is.null(one <- sap[["1"]])) cat("  One copy:", paste(one, collapse = ", "), "\n")
  if (!is.null(atm1 <- sap[["atmost1"]])) cat("  At most one copy:", paste(atm1, collapse = ", "), "\n")
  if (!is.null(zero <- sap[["0"]])) cat("  Zero copies:", paste(zero, collapse = ", "), "\n")
  cat("\n")
}