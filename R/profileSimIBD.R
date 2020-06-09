#' Simulate markers on a given IBD pattern
#'
#' This function is similar to `profileSim()` but may be used to simulate
#' genotypes for linked markers.
#'
#' @param x A `ped` object
#' @param ibdpattern A simulation output of `ibdsim`
#' @param markers A vector with names of indices of markers attached to `x`
#'
#' @return An object similar to `x`
#'
#' @examples
#' x = nuclearPed(2)
#'
#' # Simulation of IBD in the pedigree
#' s = ibdsim(x, 1, map = uniformMap(M = 1), seed = 1729)[[1]]
#' alleleSummary(s, 3:4)
#' 
#' # Attach 3 markers on chromosome 1
#' loci = list(
#'  list(afreq = c(a = 0.5, b = 0.5), chrom = 1, posCm = 20),
#'  list(afreq = c(a = 0.9, b = 0.1), chrom = 1, posCm = 30),
#'  list(afreq = c(a = 0.2, b = 0.3, c = .5), chrom = 1, posCm = 60)
#' )
#' x = setMarkers(x, loc = loci)
#'
#' # Simulate genotypes on the given IBD pattern
#' profileSimIBD(x, s)
#' 
#' 
#' @export
profileSimIBD = function(x, ibdpattern, markers = NULL) {

  if(class(ibdpattern) == "genomeSimList")
    return(lapply(ibdpattern, function(patt) profileSimIBD(x, patt, markers = markers)))
  
  if(!is.null(markers))
    x = selectMarkers(x, markers)

  a = alleleSummary(ibdpattern)
  
  nMark = nMarkers(x)
  mname = name(x, 1:nMark)
  mchr  = chrom(x, 1:nMark)
  mpos  = posCm(x, 1:nMark)
  
  if(any(is.na(mpos) | is.na(mchr)))
    stop2("All markers must have defined chromosome and position attributes")

  if(!all(mchr %in% a[, 'chrom']))
    stop2("Chromosome missing from `ibdpattern`: ", setdiff(mchr, a[, 'chrom']))
    
  # Pick rows (indices) in `a` corresponding to marker positions
  arows = vapply(seq_len(nMark), function(i) {
    chrrows = which(a[, 'chrom'] == mchr[i])
    chrrows[findInterval(mpos[i], a[chrrows, 'start'])]
  }, FUN.VALUE = 1L)
  
  # Allele columns
  acols = paste(rep(labels(x), each = 2), c('p','m'), sep = ":")
  
  # Number of founder alleles (i.e,. "different colours")
  nA = 2 * length(founders(x))
  
  # Create empty allele matrix
  alsMat.transp = matrix("0", nrow = 2*nMark, ncol = pedsize(x))
  
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
  setAlleles(x, alleles = alsMat)
}

