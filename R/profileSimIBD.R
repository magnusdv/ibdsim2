#' Simulate markers conditional on a given IBD pattern
#'
#' This function simulates genotypes for a set of markers conditional on a
#' specific underlying IBD pattern (typically produced with [ibdsim()]).
#'
#' It should be noted that the only *random* part of this function is the
#' sampling of founder alleles for each marker. Given those, all other genotypes
#' in the pedigree are determined by the underlying IBD pattern.
#'
#' @param x A `ped` object.
#' @param ibdpattern A `genomeSim()` object, typically created by [ibdsim()].
#'   (See Examples).
#' @param ids A vector of ID labels. If NULL, extracted from `ibdpattern`.
#' @param markers A vector with names or indices of markers attached to `x`.
#' @param seed An integer seed for the random number generator.
#' @param verbose A logical, by default TRUE.
#'
#' @return A copy of `x` where marker genotypes have been simulated conditional
#'   on `ibdpattern`.
#'
#' @seealso [ibdsim()], `forrel::profileSim()`.
#'
#' @examples
#' 
#' # Brother-sister pedigree
#' ped = nuclearPed(2, sex = 1:2)
#' 
#' # Alleles
#' als = letters[1:10]
#' 
#' 
#' ### Autosomal simulation
#' 
#' x = ped |> 
#'   addMarker(alleles = als, chrom = 1, posMb = 20) |> 
#'   addMarker(alleles = als, chrom = 1, posMb = 50) |> 
#'   addMarker(alleles = als, chrom = 1, posMb = 70)
#'   
#' # Simulate the underlying IBD pattern in the pedigree
#' sim = ibdsim(x, map = uniformMap(M = 1, chrom = 1), seed = 123)
#'
#' # Simulate genotypes for the sibs conditional on the given IBD pattern
#' profileSimIBD(x, sim, ids = 3:4, seed = 123)
#'
#' # With a different seed
#' profileSimIBD(x, sim, ids = 3:4, seed = 124)
#'
#'
#' ### X chromosomal simulation
#' 
#' y = ped |> 
#'   addMarker(alleles = als, chrom = "X", posMb = 1) |> 
#'   addMarker(alleles = als, chrom = "X", posMb = 50) |> 
#'   addMarker(alleles = als, chrom = "X", posMb = 100)
#'
#' simy = ibdsim(y, map = loadMap("decode19", chrom = 23), seed = 11)
#'
#' profileSimIBD(y, simy, seed = 12)
#'
#' @export
profileSimIBD = function(x, ibdpattern, ids = NULL, markers = NULL, seed = NULL, verbose = TRUE) {
  
  # Set seed if given
  if(!is.null(seed))
    set.seed(seed)
  
  if(!is.data.frame(ibdpattern) && is.list(ibdpattern))
    return(lapply(ibdpattern, function(patt) profileSimIBD(x, patt, ids = ids, markers = markers)))
  
  if(!is.null(markers))
    x = selectMarkers(x, markers)
  
  a = ibdpattern
  if(is.null(ids)) {
    ids = extractIds(a)
    if(verbose) cat("IDs extracted from provided IBD pattern:", toString(ids), "\n")
  }
  else
    a = alleleFlow(a, ids, addState = FALSE)

  if(!all(ids %in% labels(x)))
    stop2("ID label in `ibdpattern` not found in `x`: ", setdiff(ids, labels(x)))
  
  nMark = nMarkers(x)
  if(nMark == 0)
    stop2("The pedigree has no markers attached")
  
  idsInt = internalID(x, ids)
  mchr  = chrom(x, 1:nMark)
  mpos  = posMb(x, 1:nMark)
  
  if(any(is.na(mpos) | is.na(mchr)))
    stop2("All markers must have defined chromosome and position attributes")
  
  Xchrom = isXsim(a)
  
  # X-chromosome is notated as "23" in the simulation
  if(Xchrom)
    mchr[mchr == "X"] = "23"
  
  # Split a on chrom (NB: split(a, a[,'chrom']) doesn't work directly)
  aChr = lapply(split(1:nrow(a), a[, "chrom"]), function(rws) a[rws, , drop = FALSE])
  
  if(!all(mchr %in% names(aChr)))
    stop2("Chromosome missing from `ibdpattern`: ", setdiff(mchr, achr))
  
  # Allele columns
  matcols = 5 + seq_along(ids)*2L 
  patcols = matcols - 1L           
  
  # Number of founder alleles (i.e,. "different colours")
  f2 = 2 * length(founders(x))
  
  # Marker matrix template
  tmpMat = matrix(0L, nrow = pedsize(x), ncol = 2)
  
  # Fill in allele matrix one marker at a time
  for(i in seq_len(nMark)) {
    
    # IBD pattern for this marker
    achr = aChr[[mchr[i]]]
    nr = nrow(achr)
    
    interv = c(achr[, 'startMB'], achr[nr, 'endMB'])
    rw = if(nr > 1) findInterval(mpos[i], interv, all.inside = TRUE) else 1L
    ibdpat = achr[rw, patcols]
    ibdmat = achr[rw, matcols]
    
    # Marker allele frequencies
    m = x$MARKERS[[i]]
    frq = attr(m, "afreq")
    
    # Sample founder alleles
    founderAlleles = sample.int(length(frq), size = f2, replace = TRUE, prob = frq)
    
    # Ad hoc (but good enough) fix for X males
    if(Xchrom) {
      zz = ibdpat == 0
      ibdpat[zz] = ibdmat[zz]
    }
                                                     
    # Distribute alleles according to IBD pattern
    amat = tmpMat
    amat[idsInt, 1] = founderAlleles[ibdpat]
    amat[idsInt, 2] = founderAlleles[ibdmat]
    
    # Sort genotypes
    swap = amat[,1] > amat[,2]
    if(any(swap))
      amat[swap, 1:2] = amat[swap, 2:1]
    
    # Insert in marker object
    x$MARKERS[[i]][] = amat
  }
  
  x
}

