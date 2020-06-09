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
  
  if(!is.null(markers))
    x = selectMarkers(x, markers)
  
  nMark = nMarkers(x)
  mnames = name(x, 1:nMark)
  chroms = chrom(x, 1:nMark)
  cm = posCm(x, 1:nMark)
  
  if(length(unique.default(chroms)) > 1) stop2("Something wrong related to chromosomes")
  if(any(is.na(cm)))
    stop2("All markers must have defined cM positions")
  
  if(class(ibdpattern) == "genomeSimList") {
    if(length(ibdpattern) == 1)
      ibdpattern = ibdpattern[[1]]
    else
      stop2("`ibdpattern` must be a single simulation, not a list of several.")
  }
  if(length(ibdpattern) > 1)
    stop2("`ibdpattern` must be for one chromosome only, for now...")
  
  #-----------------------------
  
  a = alleleSummary(ibdpattern)
  
  # Pick rows (indices) in `a` corresponding to marker positions
  arows = findInterval(cm, a[,'start'])
  
  # Allele columns
  acols = paste(rep(labels(x), each = 2), c('p','m'), sep = ":")
  
  # Number of founder alleles (i.e,. "different colours")
  nA = 2 * length(founders(x))
  
  # Create empty allele matrix
  alsMat.transp = matrix("0", nrow = 2*nMark, ncol = pedsize(x))
  
  # Fill in allele matrix one marker at a time
  for(i in seq_len(nMark)) {
    m = mnames[i]
    freqs = afreq(x, m)
    als = names(freqs)
    frq = as.numeric(freqs)
    
    # IBD pattern for this marker
    ibdpatt = a[arows[i], acols]
    
    # Sample founder allele labels
    founderAlleles = sample(als, size = nA, replace = TRUE, prob = frq)
    
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

