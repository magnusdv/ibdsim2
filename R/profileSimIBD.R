#' Simulate markers conditional on a given IBD pattern
#'
#' This function simulates genotypes for a set of markers conditional on a
#' specific underlying IBD pattern, typically produced with [ibdsim()].
#'
#' Founder alleles are sampled independently at each marker, so linkage
#' disequilibrium is ignored. Existing genotypes and mutation models are ignored.
#' Individuals not included in `ids` have missing genotypes in the result.
#'
#' It should be noted that the only *random* part of this function is the
#' sampling of founder alleles. Given those, all other genotypes in the pedigree
#' are determined by the underlying IBD pattern and the marker positions.
#'
#' @param x A `ped` object.
#' @param ibdpattern A `genomeSim` object, or a list of such objects.
#' @param ids A vector of ID labels. If NULL, extracted from `ibdpattern`.
#' @param markers A vector with names or indices of markers attached to `x`.
#' @param seed An integer seed for the random number generator.
#' @param verbose A logical, by default TRUE.
#'
#' @return A copy of `x` with simulated genotypes, or a list of such copies if
#'   `ibdpattern` is a list.
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
profileSimIBD = function(x, ibdpattern, ids = NULL, markers = NULL, seed = NULL, 
                         verbose = TRUE) {
  
  # Set seed if given
  if(!is.null(seed))
    set.seed(seed)
  
  if(!is.null(markers))
    x = selectMarkers(x, markers)

  if(!is.data.frame(ibdpattern) && is.list(ibdpattern))
    return(lapply(ibdpattern, function(patt)
      profileSimIBD(x, patt, ids = ids, markers = NULL, verbose = verbose)))
  
  a = ibdpattern
  if(is.null(ids)) {
    ids = extractIds(a)
    if(verbose) cat("IDs extracted from provided IBD pattern:", toString(ids), "\n")
  }
  else
    a = alleleFlow(a, ids, addState = FALSE)

  if(anyNA(match(ids, labels(x))))
    stop2("ID label in `ibdpattern` not found in `x`: ", .mysetdiff(ids, labels(x)))
  
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
  
  # Allele columns
  matcols = 5 + seq_along(ids) * 2L
  patcols = matcols - 1L
  ibdcols = c(patcols, matcols)

  # Remap the founder allele codes actually used
  ibd = a[, ibdcols, drop = FALSE]
  if(anyNA(ibd) || any(ibd < 0 | ibd != as.integer(ibd)))
    stop2("Invalid founder allele code in `ibdpattern`")

  codes = sort(unique.default(as.integer(ibd[ibd > 0])))
  if(!length(codes))
    stop2("No founder allele codes found in `ibdpattern`")

  a[, ibdcols] = match(ibd, codes, nomatch = 0L)
  f2 = length(codes)

  # Split the pattern by chromosome
  aChr = lapply(split(seq_len(nrow(a)), a[, "chrom"]),
                function(rws) a[rws, , drop = FALSE])

  if(anyNA(match(mchr, names(aChr))))
    stop2("Chromosome missing from `ibdpattern`: ",
          .mysetdiff(mchr, names(aChr)))

  # Locate all markers chromosome-wise
  aRow = integer(nMark)
  for(chr in unique.default(mchr)) {
    achr = aChr[[chr]]
    idx = mchr == chr
    interv = c(achr[, "startMB"], achr[nrow(achr), "endMB"])
    aRow[idx] = findInterval(mpos[idx], interv, all.inside = TRUE)
  }
  
  # Fill in allele matrix one marker at a time
  for(i in seq_len(nMark)) {
    
    # IBD pattern for this marker
    achr = aChr[[mchr[i]]]
    rw = aRow[i]
    ibdpat = achr[rw, patcols]
    ibdmat = achr[rw, matcols]
    
    # Marker allele frequencies
    m = x$MARKERS[[i]]
    frq = attr(m, "afreq")
    
    # Sample founder alleles
    founderAlleles = if(length(frq) == 2L)
      1L + (runif(f2) > frq[1])
    else
      sample.int(length(frq), size = f2, replace = TRUE, prob = frq)
    
    # Ad hoc (but good enough) fix for X males
    if(Xchrom) {
      zz = ibdpat == 0
      ibdpat[zz] = ibdmat[zz]
    }
                                                     
    # Distribute alleles according to the IBD pattern
    m[] = 0L
    m[idsInt, 1] = founderAlleles[ibdpat]
    m[idsInt, 2] = founderAlleles[ibdmat]

    # Sort the selected genotypes
    swap = m[idsInt, 1] > m[idsInt, 2]
    if(any(swap)) {
      idx = idsInt[swap]
      m[idx, 1:2] = m[idx, 2:1]
    }

    x$MARKERS[[i]] = m
  }
  
  x
}

