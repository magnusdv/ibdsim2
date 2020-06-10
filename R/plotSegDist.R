#' Segment distribution plot
#'
#' @param segDist A data frame or a list of data frames
#' @param labels A character vector, used if `segDist` is an unnamed list of length > 1
#' @param alpha A transparency parameter for the scatter points.
#' @param ellipses A logical: Should confidence ellipses be added to the plot?
#' @param legend_inside A logical indicating wether the legend should be placed
#'   inside (default) or outside the plot window.
#' @param title,xlab,ylab Title and axis labels. Set to `NULL` to remove.
#' 
#' @examples
#' #################################################################
#' # EXAMPLE
#' # Comparison of autozygosity distributions in various individuals
#' # with the same expected inbreeding coefficient (f = 1/8)
#' #################################################################
#' 
#' G = swapSex(linearPed(2), 5)           # grandfather/granddaughter
#' G = addChildren(G, 1, 5, 1)
#' HSpat = swapSex(halfSibPed(), 5)       # paternal half sibs
#' HSpat = addChildren(HSpat, 4, 5, 1)
#' HSmat = swapSex(HSpat, 1)              # maternal half sibs
#' QHFC = quadHalfFirstCousins()          # quad half first cousins
#' QHFC = addChildren(QHFC, 9, 10, nch = 1)
#' 
#' peds = list(G = G, HSpat = HSpat, HSmat = HSmat, QHFC = QHFC)
#' plotPedList(peds, newdev = TRUE)
#' dev.off()
#' 
#' # Simulations (increase 'sims'!))
#' s = lapply(peds, ibdsim, sims = 10, map = "uniform.sex.spec")
#' 
#' # Summarise autozygous regions
#' segs = lapply(s, realisedAutozygosity, id = "leaf")
#' 
#' # Plot distributions
#' plotSegDist(segs, title = "Distribution of autozygous segments")
#' 
#' @import ggplot2
#' @importFrom ribd inbreeding
#' @export
plotSegDist = function(segDist, labels = NULL, alpha = 1, ellipses = TRUE, 
                     legend_inside = TRUE, title = NULL, 
                     xlab = "Segment count", ylab = "Average segment length (cM)") {
  
  segDistList = if(is.data.frame(segDist)) list(segDist) else segDist
  N = length(segDistList)
  
  # Labels
  if(is.null(labels))
    labels = names(segDistList)
  if(is.null(labels)) 
    labels = as.character(1:N)
  if(N > 1 && any(labels == ""))
    stop2("Only some of the elements of `segDist` are named")
  
  plotData = do.call(rbind, segDistList)
  plotData$label = factor(rep(labels, sapply(segDistList, nrow)), 
                          levels = labels)
  nLabs = nlevels(plotData$label)
  
  max.x = max(plotData$segCount)
  max.y = max(plotData$meanLength)
  
  # Create plot
  g = ggplot(data = plotData, aes_string(x = "segCount", y = "meanLength", color = "label")) + 
    geom_jitter(width = 0.25, alpha = alpha) +
    theme_bw(base_size = 15) + 
    scale_color_manual(values = ggplotColors(nLabs)) +
    labs(title = title, x = xlab, y = ylab, color = "Relationship")
  
  if(ellipses) 
    g = g + stat_ellipse(aes_string(x = "segCount"), size = 1.3)
  
  # Theoretical expectation curves
  expected = sort(unique(plotData$expected))
  genomeLen = plotData$genomeLength[1]
  
  # Remove 0 if included
  if(0 %in% expected) {
    message("No theoretical curve to draw for expectated value 0")
    expected = expected[expected > 0]
  }
  
  # Add curves to plot
  if(length(expected) > 0) {
    curveData = do.call(rbind, lapply(expected, function(v) {
      xmin = if(max.y > 0) genomeLen * v / max.y else 0
      xvec = seq(xmin, max.x, length = 50)
      data.frame(x = xvec, y = genomeLen * v / xvec, 
                 coeff = as.character(round(v, 3)))
    }))
    
    g = g + 
      geom_line(data = curveData, aes_string("x", "y", linetype = "coeff"), 
                lwd = 1, inherit.aes = FALSE) + 
      scale_linetype_manual(values = 1 + seq_along(expected), # avoid solid line
                            name = "Expected")
  }
  
  # Fix legends
  g = g + 
    theme(legend.key.width = unit(0.9, "cm")) + 
    guides(color = if(nLabs > 1) guide_legend(order = 1) else FALSE,
           linetype = guide_legend(order = 2, reverse = TRUE))
  
  if(legend_inside) 
    g = g + theme(legend.position = c(.95, .95), 
                  legend.justification = c("right", "top"))
  
  g
}

 