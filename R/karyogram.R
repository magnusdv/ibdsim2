#' Karyogram plots
#'
#' Functions for visualising IBD segments in karyograms. The `karyogram1()` and
#' `karyogram2()` functions produces karyograms illustrating the output of
#' [ibdim()] for one or two specified individuals. The actual plotting is done
#' by functions `karyoHaploid()` and `karyoDiploid()`.
#'
#' @param sim A `genomeSimList` object, typically produced by [ibdsim()].
#' @param ids A vector of one or two ID labels.
#' @param verbose A logical.
#' @param ... Further arguments passed on to `karyoHaploid()`.
#'
#' @return
#' 
#' @examples
#' x = quadHalfFirstCousins()
#' s = ibdsim(x, ids = leaves(x))
#'
#' karyogram2(s[[1]])
#' 
karyogram2 = function(sim, ids = NULL, verbose = TRUE, ...) {
  
  # IDs present in sim
  simids = extractIdsFromSegmentSummary(sim)
  
  if(is.null(ids)) {
    ids = simids
    if(verbose) 
      message("ID labels extracted from input: ", toString(ids))
  }
  
  if(length(ids) != 2)
    stop2("The input must specify exactly two ID labels: ", ids)
  
  if(!all(ids %in% simids))
    stop2("Target ID not found in input data:", setdiff(ids, simids))
  
  if(length(simids) > 2 || !"IBD" %in% colnames(sim)) {
    sim0 = segmentSummary(sim, ids = ids, addState = TRUE)
    sim = mergeSegments(sim0, by = "IBD")
  }
  
  ibd = sim[, 'IBD']
  a = sim[ibd > 0, , drop = FALSE]
  adf = as.data.frame(a)
  
  # Extract column vectors with paternal and maternal alleles
  pat1 = a[, paste(ids[1], "p", sep = ":")]
  pat2 = a[, paste(ids[2], "p", sep = ":")]
  mat1 = a[, paste(ids[1], "m", sep = ":")]
  mat2 = a[ ,paste(ids[2], "m", sep = ":")]
  
  segData = NULL
  if(any(pat1 == pat2))
    segData = rbind(segData, cbind(subset(adf, pat1 == pat2), type = 1))
  if(any(mat1 == mat2))
    segData = rbind(segData, cbind(subset(adf, mat1 == mat2), type = 2))
  if(any(pat1 == mat2))
    segData = rbind(segData, cbind(subset(adf, pat1 == mat2), type = 3))
  if(any(pat2 == mat1))
    segData = rbind(segData, cbind(subset(adf, pat2 == mat1), type = 4))
  
  segData$type = droplevels(factor(segData$type, levels = 1:4, 
                                   labels = c("Paternal", "Maternal", "Fa1 = Mo2", "Fa2 = Mo2")))
  
  karyoHaploid(segData, colorBy = "type", color = c(2,4,3,5), ...)
}


karyogram1 = function(sim, id = NULL, type = c("all", "autozygous"), verbose = TRUE, ...) {
  
  # IDs present in sim
  simids = extractIdsFromSegmentSummary(sim)
  
  if(is.null(id)) {
    id = simids
    if(verbose) 
      message("ID labels extracted from input: ", toString(id))
  }
  
  if(length(id) != 1)
    stop2("The input must specify exactly one ID label: ", id)
  
  if(!id %in% simids)
    stop2("Target ID not found in input data:", id)
  
  if(length(simids) > 2 || !"Aut" %in% colnames(sim)) {
    sim0 = segmentSummary(sim, ids = id, addState = TRUE)
    sim = mergeSegments(sim0, by = "Aut")
  }
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
#'   colouring. If NA (default), all segments will have the same colour,
#'   controlled by the `color` parameter.
#' @param color A single fill colour for all the segments, or (if `colorBy` is
#'   not NA) a named vector of colours. In the latter case, the names should
#'   include all entries in the `colorBy` column.
#' @param separate A logical; relevant only if the `colorBy` column has more
#'   than one level. If FALSE, all segments are drawn in full height. This may
#'   not be optimal if segments of different colours overlap. If TRUE the levels
#'   are drawn in separate bands on the chromosomes.
#' @param alpha A single numeric in `[0,1]` indicating colour transparency.
#' @param bgcol The background colour of the chromosomes.
#' @param title Plot title.
#'
#' @return The plot object is returned invisibly, so that additional `ggplot`
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
#' karyoHaploid(segs, color = "cyan")
#' karyoHaploid(segs, colorBy = "IBD", color = cols)
#' 
#' # Note difference if `separate = FALSE`
#' karyoHaploid(segs, colorBy = "IBD", color = cols, separate = FALSE)
#' 
#' # To see the overlap now, you can reduce alpha:
#' karyoHaploid(segs, colorBy = "IBD", color = cols, separate = FALSE, alpha = 0.7)
#'               
#' # Example showing simulated IBD segments of full siblings
#' s = ibdsim(nuclearPed(2), N = 1, ids = 3:4)
#' karyogram(s[[1]])
#' }
#'
#' 
karyoHaploid = function(segments, chrom = 1:22, colorBy = NA, color = NULL, separate = TRUE, 
                        alpha = 1, bgcol = "gray95", title = NULL, legendTitle = "IBD", 
                        prefix = "chr") {
  
  map = loadMap("decode19", chrom = chrom)
  chrlen = sapply(map, physRange)
  seqnames = paste0(prefix, sapply(map, attr, 'chrom'))
  
  genome = data.frame(chr = factor(seqnames, levels = seqnames), 
                      Mb = chrlen)
  
  segments = segments[segments[,1] %in% chrom, , drop = F]
  segments = prepare_segments(segments, colorBy)
  segments$chr = factor(paste0(prefix, segments$chr), 
                        levels = levels(genome$chr))
  
  levN = nlevels(segments$fill)
  if(is.null(color))
    color = if(!is.na(colorBy)) 1 else (1 + 1:levN)
  
  # Segments y positions
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
    ggtitle(title) +
    facet_grid(chr ~ ., switch = "y") + 
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_fill_manual(values = color) + 
    guides(fill = if(is.na(colorBy)) "none" else guide_legend(legendTitle)) +
    theme_void(base_size = 16) + 
    theme(plot.margin = margin(4, 4, 4, 4),
          plot.title = element_text(size = 16, 
                                    margin = margin(b = 10, unit = "pt")),
          strip.text.y.left = element_text(angle = 0, hjust = 1),
          legend.position = c(0.97, 0.4),
          legend.justification = c(1, 0)
          )
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
#' @param chrom The (autosomal) chromosomes to be included in the plot,
#'   given as a subset of the integers 1, 2,..., 22.
#' @param colors A vector of two colours (in any form recognisable by R). If only
#'   one colour is given it is recycled. If the vector is named, a colour legend
#'   is included in the plot, using the names as labels.
#' @param alpha A single numeric in `[0,1]` indicating colour transparency.
#' @param bgcol The background colour of the chromosomes.
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
#' karyoDiploid(pat, mat)
#' }
#'
#' 
karyoDiploid = function(paternal, maternal, chrom = 1:22, 
                         colors = c(paternal = "lightblue", maternal = "orange"),  
                         alpha = 1, bgcol = "gray95", title = NULL) {
  
  decode = loadMap("Decode", chrom = chrom)
  chrlen = sapply(decode, attr, 'length')
  seqnames = paste0("chr", chrom)
  genome = data.frame(chr = factor(seqnames, levels = seqnames), Mb = chrlen)
  
  paternal = prepare_segments(paternal, colorBy = NA)
  maternal = prepare_segments(maternal, colorBy = NA)
  
  paternal$chr = factor(paternal$chr, levels = levels(genome$chr))
  maternal$chr = factor(maternal$chr, levels = levels(genome$chr))
  
  Np = nrow(paternal)
  Nm = nrow(maternal)
   
  # Colours
  colorlabels = names(colors) %||% c("paternal", "maternal")
  colorlabels = factor(colorlabels, levels = colorlabels)
  
  paternal$fill = rep_len(colorlabels[1], Np)
  maternal$fill = rep_len(colorlabels[2], Nm)
  
  p = ggplot() + 
    theme_void(base_size = 15) +  
    ggtitle(title) +
    geom_rect(aes_string(xmin = 0, xmax = "Mb", ymin = 0, ymax = .43), 
              data = genome, fill = bgcol, col = "black") + 
    geom_rect(aes_string(xmin = 0, xmax = "Mb", ymin = 0.57, ymax = 1), 
              data = genome, fill = bgcol, col = "black") +
    facet_grid(chr ~ ., switch = "y") + 
    theme(
      plot.margin = unit(c(10, 5, 10, 5), "pt"),
      strip.text.y.left = element_text(angle = 0, margin = margin(0,0,0,0), hjust = 1),
      panel.spacing.y = unit(.25, "lines"))
  
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

# Prepare input 'segments' for plotting (extract and rename relevant columns,
# merge adjacent segments)
prepare_segments = function(segments, colorBy = NA) {
  segments = as.data.frame(segments)
  
  stopifnot(ncol(segments) >= 3, 
            is.na(colorBy) || colorBy %in% names(segments))
  
  N = nrow(segments)
  df = segments[1:3]
  names(df) = c("chr", "start", "end")
  
  # Fill colours
  df$fill = factor(if(!is.na(colorBy)) segments[[colorBy]] else rep_len(1, N))
  
  # Early return if empty
  if(N == 0) 
    return(df)
  
  # Remove 'chr' and similar
  if(grepl("chr", df$chr[1], ignore.case = TRUE))
    df$chr = sub("chr", "", sub("chrom", "", sub("chromosome", "", df$chr)))
  
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

