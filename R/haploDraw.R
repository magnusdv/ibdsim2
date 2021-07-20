#' Draw haplotypes onto a pedigree plot
#'
#' Visualise the IBD pattern of a single chromosome, by drawing haplotypes onto
#' the pedigree.
#'
#' @param x A `ped` object.
#' @param ibd A `genomeSim` object.
#' @param chrom A chromosome number, needed if `ibd` contains data from multiple
#'   chromosomes.
#' @param pos A vector recycled to the length of `labels(x)`, indicating where
#'   haplotypes should be drawn relative to the pedigree symbols: 0 = no
#'   haplotypes; 1 = below; 2 = left; 3 = above; 4 = right. By default, all are
#'   placed below.
#' @param cols A colour vector corresponding to the alleles in `ibd`.
#' @param height The haplotype height divided by the height of a pedigree
#'   symbol.
#' @param width The haplotype width divided by the width of a pedigree symbol.
#' @param sep The separation between haplotypes within a pair, given as a
#'   fraction of `width`.
#' @param dist The distance between pedigree symbols and the closest haplotype,
#'   given as a fraction of `width`.
#' @param ... Arguments passed on to `plot.ped()`.
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
#' # Restore graphics parameters
#' par(op)
#'
#' @importFrom graphics rect plot
#' @export
haploDraw = function(x, ibd, chrom = NULL, pos = 1, cols = NULL, 
                     height = 4, width = 0.5, sep = 0.75, dist = 1.5, ...) {
  
  if(!is.ped(x))
    stop2("Argument `x` must be a `ped` object")
  
  labs = labels(x)
  
  # Check that ids to be included are covered in `ibd`
  idsIBD = extractIdsFromSegmentSummary(ibd)
  if(!all(labs[pos > 0] %in% idsIBD))
    stop2("ID not found in `ibd` matrix: ", setdiff(labs[pos > 0], idsIBD))
  
  # Check for multiple chromosomes
  chrvec = ibd[, 'chrom']
  if(length(unique.default(chrvec)) > 1) {
    if(is.null(chrom))
      stop2("`ibd` contains data from multiple chromosomes. Use `chrom` to select one.")
    if(!chrom %in% chrvec)
      stop2("Unknown chromosome: ", chrom)
    ibd = ibd[chrvec == chrom, , drop = FALSE]
  }
  
  pos = rep(pos, length.out = pedsize(x))
  if(is.null(cols))
    cols = seq_len(2*length(founders(x)))
  
  # Basic pedigree plot
  p = plot(x, labs = "", keep.par = TRUE, ...)
  
  # Chromosome length
  L = sum(ibd[, 'length'])
  
  # Height/width of ped symbols
  symh = p$boxh
  symw = p$boxw

  H = height * symh
  W = width  * symw
  SEP = sep * W
  DIST = dist * W
  
  # Loop through all individuals in pedigree
  for(i in 1:pedsize(x)) {
    if(pos[i] == 0) 
      next
    
    # Center of haplo-pair
    if(pos[i] == 1) {
      X = p$x[i]
      Y = p$y[i] + symh + DIST + H/2
    }
    else if(pos[i] == 2) {
      X = p$x[i] - symw/2 - DIST - W - SEP/2
      Y = p$y[i] + symh/2
    }
    else if(pos[i] == 3) {
      X = p$x[i]
      Y = p$y[i] - DIST - H/2
    }
    else if(pos[i] == 4) {
      X = p$x[i] + symw/2 + DIST + W + SEP/2
      Y = p$y[i] + symh/2
    }
    
    id = labs[i]
    
    # Paternal haplotype
    patCol = paste0(id, ":p")
    segsPat = mergeSegments(ibd, by = patCol)
    addRect(X - SEP/2 - W/2, Y, width = W, height = H, 
            sta = segsPat[, 'start']/L, col = cols[segsPat[, patCol]])

    # Paternal haplotype
    matCol = paste0(id, ":m")
    segsMat = mergeSegments(ibd, by = matCol)
    addRect(X + SEP/2 + W/2, Y, width = W, height = H, 
            sta = segsMat[, 'start']/L, col = cols[segsMat[, matCol]])
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
