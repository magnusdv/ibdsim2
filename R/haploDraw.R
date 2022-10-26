#' Draw haplotypes onto a pedigree plot
#'
#' Visualise the IBD pattern of a single chromosome, by drawing haplotypes onto
#' the pedigree.
#'
#' @param x A `ped` object.
#' @param ibd A `genomeSim` object, typically made by [ibdsim()].
#' @param chrom A chromosome number, needed if `ibd` contains data from multiple
#'   chromosomes.
#' @param ids A vector indicating for which pedigree members haplotypes should
#'   be drawn. If NULL (default), all individuals in `ibd` are included.
#' @param pos A vector recycled to `pedsize(x)`, indicating where haplotypes
#'   should be drawn relative to the pedigree symbols: 0 = no haplotypes; 1 =
#'   below; 2 = left; 3 = above; 4 = right. By default, all are placed below.
#' @param cols A colour vector corresponding to the alleles in `ibd`.
#' @param height The haplotype height divided by the height of a pedigree
#'   symbol.
#' @param width The haplotype width, divided by the width of a pedigree symbol.
#' @param sep The separation between haplotypes within a pair, measured in
#'   pedigree symbol widths.
#' @param dist The distance between pedigree symbols and the closest haplotype,
#'   measured in pedigree symbol widths.
#' @param ... Arguments passed on to `plot.ped()`. In particular, if the
#'   haplotypes appear cropped it usually helps to increase the `margins`.
#'
#' @return None.
#'
#' @examples
#' op = par(no.readonly = TRUE)
#'
#' ###############################
#' # Example 1: A family quartet #
#' ###############################
#'
#' x = nuclearPed(2)
#' s = ibdsim(x, N = 1, map = uniformMap(M = 1), seed = 4276)
#' s[[1]]
#'
#' haploDraw(x, s[[1]], pos = c(2,4,2,4), cols = c(3,7,2,4),
#'           margins = c(2, 5, 5, 5), cex = 1.2)
#'
#'
#' ###########################
#' # Example 2: Autozygosity #
#' ###########################
#'
#' x = halfCousinPed(0, child = TRUE)
#' s = ibdsim(x, N = 1, map = uniformMap(M = 1),
#'            skipRecomb = spouses(x, 2), seed = 19499)
#' s[[1]]
#'
#' # Grey colour (8) for irrelevant founder alleles
#' haploDraw(x, s[[1]], pos = c(0,1,0,2,4,4),
#'           cols = c(8,8,3,7,8,8), margin = c(2, 2, 2, 2))
#'
#'
#' ###############################
#' # Example 3: X-chromosomal sims
#' ###############################
#'
#' x = nuclearPed(2, sex = 2:1)
#' s = ibdsim(x, N = 1, map = uniformMap(M = 1, chrom = "X"), seed = 123)
#' s[[1]]
#'
#' haploDraw(x, s[[1]], pos = c(2,4,2,4), margins = c(2, 5, 5, 5), cex = 1.2)
#'
#'
#' # Restore graphics parameters
#' par(op)
#'
#' @importFrom graphics rect plot
#' @export
haploDraw = function(x, ibd, chrom = NULL, ids = NULL, pos = 1, cols = NULL, 
                     height = 4, width = 0.75, sep = 0.75, dist = 1, ...) {
  
  if(!is.ped(x))
    stop2("Argument `x` must be a `ped` object")
  
  labs = labels(x)
  idsIBD = extractIds(ibd)
  
  if(is.null(ids))
    ids = idsIBD
  else {
    ids = as.character(ids)
    if(!all(ids %in% idsIBD))
      stop2("ID not found in `ibd` matrix: ", setdiff(ids, idsIBD))
  }
  
  if(length(pos) == 1) {
    pos = rep(pos, length(ids))
    names(pos) = ids
  }
  if(length(pos) != length(ids))
    stop2("Arguments `pos` and `ids` are incompatible")
  
  if(is.null(names(pos)))
    names(pos) = ids
  if(!all(ids %in% names(pos)))
    stop2("ID not found in `pos` vector: ", setdiff(ids, names(pos)))
  
  # Extend `pos` to all individuals
  allpos = rep(0L, pedsize(x))
  names(allpos) = labs
  allpos[ids] = pos[ids]
  
  # If `ibd` is list of 1 sim: simplify
  if(inherits(ibd, "genomeSimList") && length(ibd) == 1)
    ibd = ibd[[1]]

  Xchrom = isXsim(ibd)
  isXmale = Xchrom & (labs %in% males(x))
  
  # If `chrom` is given, check compatibility and select rows from `ibd`
  chrvec = ibd[, 'chrom']
  multipleChrom = length(unique.default(chrvec)) > 1
  if(multipleChrom && is.null(chrom))
    stop2("`ibd` contains data from multiple chromosomes. Use `chrom` to select one.")
  if(length(chrom) > 1)
    stop2("More than one chromosome selected: ", chrom)
  
  if(!is.null(chrom)) { # by now chrom has length 1
    if(chrom == 23 || chrom == "X") {
      if(!Xchrom)
        stop2("`chrom = 'X'` is indicated, but the given simulation is not X-chromosomal")
      chrom = 23
    }
    if(!chrom %in% chrvec)
      stop2("Chromosome not present in the simulation: ", chrom)
    ibd = ibd[chrvec == chrom, , drop = FALSE]
  }
  
  if(is.null(cols))
    cols = seq_len(2*length(founders(x)))
  
  # Basic pedigree plot
  p = plot(x, labs = "", keep.par = TRUE, ...)
  
  # Chromosome length
  L = sum(ibd[, 'endMB'] - ibd[, 'startMB'])
  
  # Height/width of ped symbols
  symh = p$boxh
  symw = p$boxw

  H = height * symh
  W = width  * symw
  SEP = sep * symw
  DIST = dist * symw
  
  # Loop through all individuals in pedigree
  for(i in 1:pedsize(x)) {
    if(allpos[i] == 0) 
      next
    
    # Center of haplo-pair
    if(allpos[i] == 1) {
      X = p$x[i]
      Y = p$y[i] + symh + DIST + H/2
    }
    else if(allpos[i] == 2) {
      X = p$x[i] - symw/2 - DIST - W - SEP/2
      Y = p$y[i] + symh/2
    }
    else if(allpos[i] == 3) {
      X = p$x[i]
      Y = p$y[i] - DIST - H/2
    }
    else if(allpos[i] == 4) {
      X = p$x[i] + symw/2 + DIST + W + SEP/2
      Y = p$y[i] + symh/2
    }
    
    id = labs[i]
    
    if(isXmale[i]) { # draw maternal haplotype only
      matCol = paste0(id, ":m")
      segsMat = mergeSegments(ibd, by = matCol)
      addRect(X, Y, width = W, height = H, sta = segsMat[, 'startMB']/L, col = cols[segsMat[, matCol]])
    }
    else {
      # Paternal haplotype
      patCol = paste0(id, ":p")
      segsPat = mergeSegments(ibd, by = patCol)
      addRect(X - SEP/2 - W/2, Y, width = W, height = H, 
              sta = segsPat[, 'startMB']/L, col = cols[segsPat[, patCol]])
  
      # Maternal haplotype
      matCol = paste0(id, ":m")
      segsMat = mergeSegments(ibd, by = matCol)
      addRect(X + SEP/2 + W/2, Y, width = W, height = H, 
              sta = segsMat[, 'startMB']/L, col = cols[segsMat[, matCol]])
    }
  }
}


# Function for drawing a single haplotype
addRect = function(xmid, ymid, width, height, sta = 0, col = seq_along(sta)+1) {
  stopifnot(length(col) == length(sta))
  bottom = ymid - height/2
  sto = c(sta[-1], 1)
  for(i in seq_along(sta))
    rect(xleft = xmid - width/2, ybottom = bottom + height*sta[i],
         xright = xmid + width/2, ytop = bottom + height*sto[i], 
         col = col[i])
}
