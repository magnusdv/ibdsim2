#' Simulate markers on a given IBD pattern
#'
#' This function simulates genotypes for a set of markers, conditional on a
#' specific underlying IBD pattern.
#'
#' It should be noted that the only *random* part of this function is the
#' selection of founder alleles for each marker. Given those, all other
#' genotypes in the pedigree are determined by the underlying IBD pattern.
#'
#' @param x A `ped` object.
#' @param ibdpattern A `genomeSim()` object, typically created by [ibdsim()].
#'   (See Examples).
#' @param ids A vector of ID labels. If NULL, all members of `x` are included.
#' @param markers A vector with names or indices of markers attached to `x`.
#' @param seed An integer seed for the random number generator.
#'
#' @return An object similar to `x`. but with simulated genotypes.
#'
#' @seealso [ibdsim()]
#'
#' @examples
#' # A pedigree with two siblings
#' x = nuclearPed(2)
#'
#' # Attach 3 linked markers on chromosome 1
#' pos = c(20, 50, 70)   # marker positions in megabases
#' mlist = lapply(pos, function(i)
#'   marker(x, alleles = letters[1:10], chrom = 1, posMb = i))
#' x = setMarkers(x, mlist)
#'
#' # Simulate the underlying IBD pattern in the pedigree
#' s = ibdsim(x, 1, map = uniformMap(M = 1, chrom = 1), seed = 123)[[1]]
#'
#' # Simulate genotypes for the sibs conditional on the given IBD pattern
#' profileSimIBD(x, s, ids = 3:4, seed = 123)
#'
#' # With a different seed
#' profileSimIBD(x, s, ids = 3:4, seed = 124)
#'
#' @export
profileSimIBD = function(x, ibdpattern, ids = NULL, markers = NULL, seed = NULL) {

  if(!is.data.frame(ibdpattern) && is.list(ibdpattern))
    return(lapply(ibdpattern, function(patt) profileSimIBD(x, patt, markers = markers)))
  
  if(!is.null(markers))
    x = selectMarkers(x, markers)

  a = ibdpattern
  if(is.null(ids))
    ids = extractIds(a)
  else
    a = alleleFlow(a, ids, addState = FALSE)
  
  if(!all(ids %in% labels(x)))
    stop2("ID label in `ibdpattern` not found in `x`: ", setdiff(ids, labels(x)))
  
  achr = a[, 'chrom']
  
  nMark = nMarkers(x)
  mchr  = chrom(x, 1:nMark)
  mpos  = posMb(x, 1:nMark)
  
  if(any(is.na(mpos) | is.na(mchr)))
    stop2("All markers must have defined chromosome and position attributes")

  if(!all(mchr %in% achr))
    stop2("Chromosome missing from `ibdpattern`: ", setdiff(mchr, achr))
    
  # Pick rows (indices) in `a` corresponding to marker positions
  arows = vapply(seq_len(nMark), function(i) {
    chrrows = which(achr == mchr[i])
    if(length(chrrows) == 1)
      chrrows
    else
      chrrows[findInterval(mpos[i], a[chrrows, 'start'], all.inside = TRUE)]
  }, FUN.VALUE = 1L)
  
  # Allele columns
  acols = 4 + seq_len(2 * length(ids))
  
  # Number of founder alleles (i.e,. "different colours")
  nA = 2 * length(founders(x))
  
  # Create empty allele matrix
  alsMat.transp = matrix("0", nrow = 2*nMark, ncol = length(ids))
  
  # Set seed if given
  if(!is.null(seed))
    set.seed(seed)
  
  # Fill in allele matrix one marker at a time
  for(i in seq_len(nMark)) {
    # IBD pattern for this marker
    ibdpatt = a[arows[i], acols]

    # Sample founder allele labels
    frqvec = afreq(x, i)
    als = names(frqvec)
    founderAlleles = sample(als, size = nA, replace = TRUE, prob = frqvec)

    # Distribute alleles according to IBD pattern
    allAlleles = founderAlleles[ibdpatt]
    
    # Insert in transposed matrix
    alsMat.transp[c(2*i-1, 2*i), ] = allAlleles
  }
  
  # Transpose back
  alsMat = t.default(alsMat.transp)
  
  # Attach and return
  y = setAlleles(x, ids = ids, alleles = alsMat)
  y = sortGenotypes(y)
  
  y
}

