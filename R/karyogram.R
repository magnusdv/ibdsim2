#' Karyogram plots
#'
#' Functions for visualising IBD segments in karyograms. The `karyogram1()` and
#' `karyogram2()` functions produces karyograms illustrating the output of
#' [ibdsim()] for one or two specified individuals. The actual plotting is done
#' by functions `karyoHaploid()` and `karyoDiploid()`.
#'
#' @param sim A `genomeSim` object, or a `genomeSimList` of length 1, e.g.
#'   produced by `ibdsim(..., N = 1)`.
#' @param ids A vector of one or two ID labels.
#' @param verbose A logical.
#' @param ... Further arguments passed on to `karyoHaploid()`.
#'
#' @return A plot object returned invisibly.
#'
#' @examples
#' \donttest{
#' x = quadHalfFirstCousins()
#' s = ibdsim(x, seed = 1729)
#' # karyogram2(s, ids = leaves(x), title = "QHFC")
#' }
#' 
karyogram2 = function(sim, ids = NULL, verbose = TRUE, ...) {
  
  if(inherits(sim, "genomeSimList") && length(sim) == 1)
    sim = sim[[1]]
  
  # IDs present in sim
  simids = extractIds(sim)
  
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
    sim0 = alleleFlow(sim, ids = ids, addState = TRUE)
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
  
  karyoHaploid(segData, colBy = "type", col = c(2,4,3,5), ...)
}


karyogram1 = function(sim, id = NULL, type = c("all", "autozygous"), verbose = TRUE, ...) {
  
  # IDs present in sim
  simids = extractIds(sim)
  
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
    sim0 = alleleFlow(sim, ids = id, addState = TRUE)
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
#'   possibly a column indicated by `colBy`.
#' @param chrom A vector indicating which chromosomes to include. 
#' @param colBy A character vector naming the columns to be used for
#'   colouring. If NULL (default), all segments have the same colour.
#' @param col A single fill colour for all the segments, or (if `colBy` is
#'   used) a named vector of colours. In the latter case, the names should
#'   include all entries in the `colBy` column.
#' @param separate A logical; relevant only if the `colBy` column has more
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
#'                   IBD = c("paternal","maternal"))
#' cols = c(paternal = "blue", maternal = "red")
#'
#' karyoHaploid(segs, col = "cyan")
#' karyoHaploid(segs, colBy = "IBD", col = cols)
#'
#' # Note difference if `separate = FALSE`
#' karyoHaploid(segs, colBy = "IBD", col = cols, separate = FALSE)
#'
#' # Reduce alpha to see the overlaps:
#' karyoHaploid(segs, colBy = "IBD", col = cols, separate = FALSE, alpha = 0.7)
#'
#' }
#'
#' @export
karyoHaploid = function(segments, chrom = 1:22, colBy = NULL, col = NULL, separate = TRUE, 
                        alpha = 1, bgcol = "gray92", title = NULL) {
  
  # Genome map
  map = loadMap("decode19", chrom = chrom)
  seqnames = sapply(map, attr, 'chrom')
  genome = data.frame(chrom = factor(seqnames, levels = seqnames), 
                      Mb = sapply(map, attr, "physEnd"))
  
  segments = prepare_segments(segments, chrom = chrom, colBy = colBy)
  segments$chrom = factor(segments$chrom, levels = levels(genome$chrom))
  
  if(is.null(col))
    col = if(is.null(colBy)) 2 else (1 + seq_along(levels(segments$fill)))

  # Segments y positions
  if(separate) {
    heig = 1/segments$gsize
    segments$ymin = (segments$gidx - 1) * heig
    segments$ymax = segments$gidx * heig
  }
  else {
    segments$ymin = 0
    segments$ymax = 1
  }
  
  # Build plot object
  ggplot() + 
    geom_rect(aes_(xmin = 0, xmax = ~Mb, ymin = 0, ymax = 1), 
              data = genome, fill = bgcol, col = "black") + 
    geom_rect(aes_(xmin = ~start, xmax = ~end, ymin = ~ymin, ymax = ~ymax, 
                   fill = ~fill), data = segments, color = 1, alpha = alpha) +
    ggtitle(title) +
    facet_grid(chrom ~ ., switch = "y") + 
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_fill_manual(values = col) + 
    guides(fill = guide_legend(colBy[1])) +
    theme_void(base_size = 16) + 
    theme(plot.margin = margin(4, 4, 4, 4),
          plot.title = element_text(size = 16, 
                                    margin = margin(b = 10, unit = "pt")),
          strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5),
          legend.position = c(0.97, 0),
          legend.justification = c(1, 0),
          panel.spacing.y = unit(0.3, "lines"))
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
#' @param col A vector of two colours (in any form recognisable by R). If only
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
                         col = c(paternal = "lightblue", maternal = "orange"),  
                         alpha = 1, bgcol = "gray95", title = NULL) {
  
  map = loadMap("decode19", chrom = chrom)
  chrlen = sapply(map, attr, "physEnd")
  seqnames = paste0("chr", chrom)
  genome = data.frame(chrom = factor(seqnames, levels = seqnames), Mb = chrlen)
  
  paternal = prepare_segments(paternal, chrom = chrom, colBy = NA)
  maternal = prepare_segments(maternal, chrom = chrom, colBy = NA)
  
  paternal$chrom = factor(paternal$chrom, levels = levels(genome$chrom))
  maternal$chrom = factor(maternal$chrom, levels = levels(genome$chrom))
  
  Np = nrow(paternal)
  Nm = nrow(maternal)
   
  # Colours
  colourlabels = names(col) %||% c("paternal", "maternal")
  colourlabels = factor(colourlabels, levels = colourlabels)
  
  paternal$fill = rep_len(colourlabels[1], Np)
  maternal$fill = rep_len(colourlabels[2], Nm)
  
  p = ggplot() + 
    theme_void(base_size = 15) +  
    ggtitle(title) +
    geom_rect(aes_string(xmin = 0, xmax = "Mb", ymin = 0, ymax = .43), 
              data = genome, fill = bgcol, col = "black") + 
    geom_rect(aes_string(xmin = 0, xmax = "Mb", ymin = 0.57, ymax = 1), 
              data = genome, fill = bgcol, col = "black") +
    facet_grid(chrom ~ ., switch = "y") + 
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
    scale_fill_manual(values = col, drop = FALSE) +
    labs(fill = NULL)
      
  p
}

# Prepare input 'segments' for plotting (extract and rename relevant columns,
# merge adjacent segments)
prepare_segments = function(segments, chrom = 1:22, colBy = NULL) {
  segments = as.data.frame(segments)
  
  if(ncol(segments) < 3)
    stop2("Data frame `segments` must have at least 3 columns (chrom, start, end)")
  if(!is.null(colBy) && !all(colBy %in% names(segments)))
    stop2("Unknown column name: ", setdiff(colBy, names(segments)))
  
  # Numeric chrs
  chrNo = segments[,1] = sub("chr", "", sub("chrom", "", sub("chromosome", "", segments[,1])))
  
  df = segments[chrNo %in% chrom, , drop = FALSE]
  N = nrow(df)
  
  # Fill colours
  if(is.null(colBy))
    fill = rep_len(1, N)
  else
    fill = apply(df[colBy], 1, paste, collapse = "-")
  
  df$fill = factor(fill)
  
  names(df)[1:3] = c("chrom", "start", "end")
  df = df[c("chrom", "start", "end", "fill")]
  
  # Early return if empty
  if(N == 0) 
    return(df)
  
  # Sort
  if(N > 1)
    df = df[order(as.numeric(df$chrom), df$start, df$end), , drop = FALSE]
  
  if(max(df$end) > 250e3) {
    df$start = df$start/1e6
    df$end = df$end/1e6
    message("Converting positions to Mb by diving by 1e6")
  } else if(max(df$end) > 250) {
    df$start = df$start/1e3
    df$end = df$end/1e3
    message("Converting positions to Mb by diving by 1000")
  }
  
  # Merge overlapping segments with the same colour
  df = mergeSegments(df, by = c("chrom", "fill"), checkAdjacency = TRUE)
  rownames(df) = NULL
  
  groupOverlaps(df)
}

