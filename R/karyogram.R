
# Prepare input 'segments' for plotting (extract and rename relevant columns,
# merge adjacent segments)
prepare_segments = function(segments, colorBy = NA) {
  segments = as.data.frame(segments)
  
  stopifnot(ncol(segments) >= 3, 
            is.na(colorBy) || colorBy %in% names(segments))
  
  N = nrow(segments)
  df = segments[1:3]
  names(df) = c("chr", "start", "end")
  
  # Fill colors
  df$fill = factor(if(!is.na(colorBy)) segments[[colorBy]] else rep_len(1, N))
  
  # Early return if empty
  if(N == 0) 
    return(df)
  
  # Add 'chr' if necessary
  if(df$chr[1] %in% 1:22) 
    df$chr = paste0("chr", df$chr)
  
  if(max(df$end) > 250e3) {
    df$start = df$start/1e6
    df$end = df$end/1e6
    message("Converting positions to Mb by diving by 1e6")
  }
  else if(max(df$end) > 250) {
    df$start = df$start/1e3
    df$end = df$end/1e3
    message("Converting positions to Mb by diving by 1000")
  }
  
  # Merge overlapping segments with the same color
  df = mergeConsecutiveSegments(df, mergeBy = c("chr", "fill"), segLength = "length")
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
#'   possibly a column indicated by `colorBy`.
#' @param colorBy The name of a single column of `segments`, to be used for
#'   coloring. If NA (default), all segments will have the same color,
#'   controlled by the `color` parameter.
#' @param color A single fill color for all the segments, or (if `colorBy` is
#'   not NA) a named vector of colors. In the latter case, the names should
#'   include all entries in the `colorBy` column.
#' @param separate A logical; relevant only if the `colorBy` colomn has more
#'   than one level. If FALSE, all segments are drawn in full height. This may
#'   not be optimal if segments of different colors overlap. If TRUE the levels
#'   are drawn in separate bands on the chromosomes.
#' @param alpha A single numeric in `[0,1]` indicating color transparency.
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
#' segs = data.frame(chrom = c(1,4,5,5,10,10),
#'                   start = c(100,50,20,80,10,50),
#'                   end = c(120,100,25,100,70,120),
#'                   IBD = c("cousin1","cousin2"))
#' cols = c(cousin1 = "blue", cousin2 = "red")
#' 
#' karyo_haploid(segs, color = "cyan")
#' karyo_haploid(segs, colorBy = "IBD", color = cols)
#' 
#' # Note difference if `separate = FALSE`
#' karyo_haploid(segs, colorBy = "IBD", color = cols, separate = FALSE)
#' 
#' # To see the overlap now, you can reduce alpha:
#' karyo_haploid(segs, colorBy = "IBD", color = cols, separate = FALSE, alpha = 0.7)
#'               
#' # Example showing simulated IBD segments of full siblings
#' x = nuclearPed(2)
#' s = ibdsim(x, sims = 1)[[1]]
#' a = as.data.frame(alleleSummary(s, 3:4))
#' a$status = "No IBD"
#' a$status[a$IBD == 1 & a$`3:p` == a$`4:p`] = "Paternal"
#' a$status[a$IBD == 1 & a$`3:m` == a$`4:m`] = "Maternal"
#' a$status[a$IBD == 2] = "Pat & mat"
#' a$status = as.factor(a$status)
#'
#' karyo_haploid(a, colorBy = "status", separate = FALSE)
#' }
#'
#' @export
karyo_haploid = function(segments, colorBy = NA, color = "black", separate = TRUE, 
                         alpha = 1, bgcol = "gray98", title = NULL) {
  
  decode = loadMap("Decode", chrom = 1:22)
  chrlen = sapply(decode, attr, 'length')
  seqnames = paste0("chr", 1:22)
  genome = data.frame(chr = factor(seqnames, levels = seqnames), 
                      Mb = chrlen)
  
  segments = prepare_segments(segments, colorBy)
  segments$chr = factor(segments$chr, levels = levels(genome$chr))
  
  levN = nlevels(segments$fill)
  if(!is.na(colorBy) && length(color) != levN) 
    color = 1:levN
  
  # segments y positions
  if(separate) {
    heig = levN + 1 - as.integer(segments$fill)
    segments$ymin = (heig - 1)/levN
    segments$ymax = heig/levN
  }
  else {
    segments$ymin = 0
    segments$ymax = 1
  }
  
  # Build plot object
  p = ggplot() + 
    geom_rect(data = genome, 
              aes_string(xmin = 0, xmax = "Mb", ymin = 0, ymax = 1), 
              fill = bgcol, col = "black") + 
    geom_rect(data = segments, 
              aes_string(xmin = "start", xmax = "end", ymin = "ymin", ymax = "ymax", 
                         fill = "fill"), 
              col = "black", alpha = alpha) +
    facet_grid(chr ~ ., switch = "y") + 
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_fill_manual(values = color) + 
    theme_void() + 
    theme(plot.margin = margin(4, 4, 4, 4),
          plot.title = element_text(size = 16, 
                                    margin = margin(b = 10, unit = "pt")),
          strip.text.y = element_text(angle = 180, hjust = 1),
          legend.position = c(1, 0.4),
          legend.justification = c(1, 0),
          legend.text = element_text(size = 14)) +
    labs(fill = NULL, title = title)
  
  if(is.na(colorBy)) 
    p = p + guides(fill = "none")
  
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
#' @param chromosomes The (autosomal) chromosomes to be included in the plot,
#'   given as a subset of the integers 1, 2,..., 22.
#' @param colors A vector of two colors (in any form recognisable by R). If only
#'   one color is given it is recycled. If the vector is named, a color legend
#'   is included in the plot, using the names as labels.
#' @param alpha A single numeric in `[0,1]` indicating color transparency.
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
#' pat = data.frame(chrom = c(1,4,5,5,10,10), start = c(100,50,20,80,10,80),
#'                  end = c(120,100,25,100,70,120))
#' mat = data.frame(chrom = c(2,4,5,5,10), start = c(80,50,10,80,50),
#'                  end = c(120,100,35,100,120))
#' karyo_diploid(pat, mat)
#' }
#'
#' @export
karyo_diploid = function(paternal, maternal, chromosomes = 1:22, 
                         colors = c(paternal = "lightblue", maternal = "orange"),  
                         alpha = 1, bgcol = "gray99", title = NULL) {
  
  decode = loadMap("Decode", chrom = chromosomes)
  chrlen = sapply(decode, attr, 'length')
  seqnames = paste0("chr", chromosomes)
  genome = data.frame(chr = factor(seqnames, levels = seqnames), Mb = chrlen)
  
  paternal = prepare_segments(paternal, colorBy = NA)
  maternal = prepare_segments(maternal, colorBy = NA)
  
  paternal$chr = factor(paternal$chr, levels = levels(genome$chr))
  maternal$chr = factor(maternal$chr, levels = levels(genome$chr))
  
  Np = nrow(paternal)
  Nm = nrow(maternal)
   
  # colors
  colorlabels = names(colors) %||% c("paternal", "maternal")
  colorlabels = factor(colorlabels, levels = colorlabels)
  
  paternal$fill = rep_len(colorlabels[1], Np)
  maternal$fill = rep_len(colorlabels[2], Nm)
  
  p = ggplot() + 
    theme_void(base_size = 15) + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +  
    ggtitle(title) +
    geom_rect(aes_string(xmin = 0, xmax = "Mb", ymin = 0, ymax = .43), 
              data = genome, fill = bgcol, col = "black") + 
    geom_rect(aes_string(xmin = 0, xmax = "Mb", ymin = 0.57, ymax = 1), 
              data = genome, fill = bgcol, col = "black") +
    facet_grid(chr ~ ., switch = "y") + 
    theme(strip.text.y = element_text(angle = 180),
          panel.spacing.y = unit(.2, "lines"))
  
  if(Np > 0)
    p = p + geom_rect(data = paternal, 
                      aes_string(xmin = "start", xmax = "end", 
                                 ymin = 0.57, ymax = 1, fill = "fill"), 
                      col = "black", alpha = alpha)
  if(Nm > 0)
    p = p +  geom_rect(data = maternal, 
                       aes_string(xmin = "start", xmax = "end", 
                                  ymin = 0, ymax = .43, fill = "fill"), 
                       col = "black", alpha = alpha)
  
  p = p + 
    scale_fill_manual(values = colors, drop = FALSE) +
    labs(fill = NULL)
      
  p
}

