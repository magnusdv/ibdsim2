#' Allele sharing summary
#'
#' This function facilitates downstream analysis of simulations produced by
#' [ibdsim()]. It summarises a single genome simulation by describing the allele
#' flow through the pedigree.
#'
#' @param x An object of class `genomeSim`, i.e. a list of simulated
#'   chromosomes. Each chromosome is a list, with one entry for each individual.
#'   Each of these entries is a list of two matrices (one for each strand). The
#'   matrices have 2 columns (start position; allele) and one row for each
#'   segment unbroken by recombination.
#' @param ids A vector of ID labels. If missing, all individuals are included.
#'
#' @return A numerical matrix. Each row corresponds to a chromosomal segment.
#'   The first 4 columns describe the segment (chromosome, start, end, length),
#'   and are followed by two columns (paternal allele, maternal allele) for each
#'   of the selected individuals. If `length(ids) == 2` two additional columns
#'   are added:
#'
#'   * `IBD` : The IBD status of each segment (= number of alleles shared
#'   identical by descent). For a given segment, the IBD status is either 0, 1,
#'   2 or NA. If either individual is inbred, they may be autozygous in a
#'   segment, in which case the IBD status is reported as NA. With inbred
#'   individuals the `Sigma` column (see below) is more informative than the
#'   `IBD` column.
#'
#'   * `Sigma` : The condensed identity ("Jacquard") state of each segment,
#'   given as an integer in the range 1-9. The numbers correspond to the
#'   standard ordering of the condensed states. In particular, for non-inbred
#'   individuals the states 9, 8, 7 correspond to IBD status 0, 1, 2
#'   respectively.
#'
#' @examples
#' ### Sibling simulation (3 sims of chromosomes 1 and 2)
#' x = pedtools::nuclearPed(2)
#' sim = ibdsim(x, sims = 3, chromosomes = 1:2)
#' 
#' sim1 = sim[[1]] # the first simulation
#'
#' alleleSummary(sim1) # Summary of all individuals
#' alleleSummary(sim1, ids = 3:4) # Summary of the siblings
#' alleleSummary(sim1, ids = 1) # Summary of father. Trivial!
#' 
#' # Full sib mating: all 9 states are possible
#' y = pedtools::fullSibMating(1)
#' sim = ibdsim(y, sims = 1, chrom = 1, seed = 36)[[1]]
#' a = alleleSummary(sim, ids = 5:6)
#' 
#' stopifnot(setequal(a[, 'Sigma'], 1:9))
#' 
#' @importFrom pedtools internalID
#' @export
alleleSummary = function(x, ids) {
  if(!inherits(x, "genomeSim")) 
    stop2("The first argument is not a `genomeSim` object")
  
  ped = attr(x, "pedigree")
  if (missing(ids)) 
    ids = labels(ped)
  
  allele.colnames = paste0(rep(ids, each = 2), c(":p", ":m"))

  ids_int = internalID(ped, ids) 
  
  each.chrom = lapply(x, function(y) {
    haplos = unlist(y[ids_int], recursive = FALSE)
    breaks = unlist(lapply(haplos, function(m) m[-1, 1]))
    breaks = c(0, .sortDouble(breaks[!duplicated(breaks)]))
    alleles = vapply(haplos, pos2allele, posvec = breaks, FUN.VALUE = breaks)
    if (length(breaks) == 1) # Ad hoc fix, since vapply simplifies if FUN.VALUE has length 1
      dim(alleles) = c(1, length(allele.colnames))
    colnames(alleles) = allele.colnames

    stops = c(breaks[-1], attr(y, "length_Mb"))
    chrom = rep.int(attr(y, "chromosome"), length(breaks))
    cbind(chrom = chrom, start = breaks, end = stops, length = stops - breaks, alleles)
  })
  
  res = do.call(rbind, each.chrom)
  
  # Extra columns when "ids" is a pair
  if(length(ids) == 2) {
    # Jacquard state
    Sigma = jacquardState(a1 = res[, 5], a2 = res[, 6], b1 = res[, 7], b2 = res[, 8])
    
    # IBD state (0 <-> 9; 1 <-> 8; 2 <-> 7; otherwise NA)
    IBD = 9 - Sigma
    IBD[IBD > 2] = NA
    
    res = cbind(res, IBD = IBD, Sigma = Sigma)
  }
  
  res
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

