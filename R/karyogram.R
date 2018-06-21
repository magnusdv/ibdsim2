
# Prepare input 'segments' for plotting (extract and rename relevant columns,
# merge adjacent segments)
prepare_segments = function(segments, colorBy=NA) {
  segments = as.data.frame(segments)
  df = segments[1:3]
  names(df) = c("chr", "start", "end")
  
  # fill colors
  if(!is.na(colorBy)) {
    stopifnot(colorBy %in% names(segments))
    df$fill = as.factor(segments[[colorBy]])
  }
  else df$fill = factor(1)
  
  # merge overlapping segments with the same color
  N = nrow(df)
  if(N > 1) {
    join = df$end[-N] == df$start[-1] & df$fill[-N] == df$fill[-1] & df$chr[-N] == df$chr[-1]
    for(i in which(join)) 
      df$start[i+1] = df$start[i]
    df = df[!join, ]
  }
  
  # add 'chr' if neccessary
  if(df$chr[1] %in% 1:22) 
    df$chr = paste0("chr", df$chr)
  df
}


#' Haploid karyogram
#'
#' Show chromosomal segments in a haploid karyogram
#'
#' @param segments A data.frame (or an object coercible to data.frame)
#'   containing the segments to be shown on the karyogram. The first three
#'   columns must contain chromosome (with or without "chr" prefix), start
#'   position and stop position (in Mb). Any further columns are ignored, except
#'   possibly a column indicated by \code{colorBy}.
#' @param colorBy The name of a single column of \code{segments}, to be used for
#'   coloring. If NA (default), all segments will have the same color,
#'   controlled by the \code{color} parameter.
#' @param color A single fill color for all the segments, or (if \code{colorBy}
#'   is not NA) a named vector of colors. In the latter case, the names should
#'   include all entries in the \code{colorBy} column.
#' @param alpha A single numeric in [0,1] indicating color transparency.
#' @param bgcol The background color of the chromosomes.
#' @param title Plot title.
#'
#' @return The plot object is returned invisibly, so that additional ggplot
#'   layers may be added if needed.
#' @import ggplot2
#'
#' @examples
#' 
#' \dontrun{
#' segs = data.frame(chrom = c(1,4,5,5,10,10), start=c(100,50,20,80,10,50), 
#'                   end = c(120,100,25,100,70,120), IBD=c("cousin1","cousin2"))
#' karyo_haploid(segs, color="cyan")
#' karyo_haploid(segs, colorBy="IBD", color=c(cousin1="blue", cousin2="red"))
#' 
#' # To see overlaps, reduce alpha:
#' karyo_haploid(segs, colorBy="IBD", color=c(cousin1="blue", cousin2="red"), alpha=0.6)
#'
#' # Example showing simulated IBD segments of full siblings
#' x = pedtools::nuclearPed(2)
#' s = ibdsim(x, sims=1)
#' a = tibble::as_tibble(alleleSummary(s[[1]], 3:4, ibd.status=TRUE))
#' karyo_haploid(subset(a,ibd>0), colorBy="ibd", color=c("1"="blue", "2"="red"), alpha=.3)
#' }
#' 
#' @export
karyo_haploid = function(segments, colorBy=NA, color="black", alpha=1, bgcol="gray99", title=NULL) {
  
  decode = loadMap("Decode", chrom=1:22)
  chrlen = sapply(decode, attr, 'length')
  seqnames = paste0("chr", 1:22)
  genome = data.frame(chr=factor(seqnames, levels=seqnames), Mb=chrlen)
  
  segments = prepare_segments(segments, colorBy)
  segments$chr = factor(segments$chr, levels=levels(genome$chr))
  if(!is.na(colorBy) && length(color) != nlevels(segments$fill)) 
    color = 1:nlevels(segments$fill)
  
  p = ggplot() + theme_void() + theme(strip.text.y = element_text(angle = 180)) +
    geom_rect(data = genome, aes_string(xmin=0, xmax="Mb", ymin=0, ymax=1), fill=bgcol, col="black") + 
    geom_rect(data = segments, aes_string(xmin="start", xmax="end", ymin=0, ymax=1, fill="fill"), col="black", alpha=alpha) +
    facet_grid(chr~., switch="y") + labs(fill=NULL) +labs(caption="Simulation by ibdsim2") +
    scale_fill_manual(values = color)
  if(is.na(colorBy)) 
    p = p + guides(fill="none")
  if(!is.null(title)) 
    p = p + ggtitle(title)
  p
}


#' Diploid karyogram
#'
#' Show chromosomal segments in a diploid karyogram
#'
#' @param paternal,maternal data.frames (or objects coercible to data.frames)
#'   containing the segments to be shown on the paternal and maternal strands of
#'   the karyogram. The first three columns must contain chromosome (with or
#'   without "chr" prefix), start position and stop position (in Mb). Column
#'   names are ignored, as well as any further columns.
#' @param colors A vector of two colors (in any form recognisable by R). If only
#'   one color is given it is recycled. If the vector is named, a color legend
#'   is included in the plot, using the names as labels.
#' @param alpha A single numeric in [0,1] indicating color transparency.
#' @param bgcol The background color of the chromosomes.
#' @param title Plot title.
#'
#' @return The plot object is returned invisibly, so that additional ggplot
#'   layers may be added if needed.
#' @import ggplot2
#'
#' @examples
#' 
#' \dontrun{
#' pat = data.frame(chrom = c(1,4,5,5,10,10), start=c(100,50,20,80,10,80),
#'                  end = c(120,100,25,100,70,120))
#' mat = data.frame(chrom = c(2,4,5,5,10), start=c(80,50,10,80,50),
#'                  end = c(120,100,35,100,120))
#' karyo_diploid(pat, mat, color=c(paternal="lightblue", maternal="orange"), alpha=0.8)
#' }
#'
#' @export
karyo_diploid = function(paternal, maternal, colors=NA, alpha=1, bgcol="gray99", title=NULL) {
  
  decode = loadMap("Decode", chrom=1:22)
  chrlen = sapply(decode, attr, 'length')
  seqnames = paste0("chr", 1:22)
  genome = data.frame(chr=factor(seqnames, levels=seqnames), Mb=chrlen)
  
  paternal = prepare_segments(paternal, colorBy=NA)
  maternal = prepare_segments(maternal, colorBy=NA)
  
  paternal$chr = factor(paternal$chr, levels=levels(genome$chr))
  maternal$chr = factor(maternal$chr, levels=levels(genome$chr))
  
  # colors
  colors = rep(colors, length=2)
  labels = names(colors)
  names(colors) = 1:2
  paternal$fill = factor(1, levels=1:2)
  maternal$fill = factor(2, levels=1:2)
  
  p = ggplot() + theme_void() + theme(strip.text.y = element_text(angle = 180)) +
    theme(panel.spacing.y = unit(.2, "lines")) +
    geom_rect(data = genome, aes_string(xmin=0, xmax="Mb", ymin=0, ymax=.43), fill=bgcol, col="black") + 
    geom_rect(data = genome, aes_string(xmin=0, xmax="Mb", ymin=0.57, ymax=1), fill=bgcol, col="black") + 
    geom_rect(data = paternal, aes_string(xmin="start", xmax="end", ymin=.57, ymax=1, fill="fill"), col="black", alpha=alpha) +
    geom_rect(data = maternal, aes_string(xmin="start", xmax="end", ymin=0, ymax=.43, fill="fill"), col="black", alpha=alpha) +
    facet_grid(chr~., switch="y") + labs(fill=NULL) +labs(caption="Simulation by ibdsim2") +
    scale_fill_manual(values = colors, labels=labels)
  if(!is.null(title)) 
    p = p + ggtitle(title)
  p
}

