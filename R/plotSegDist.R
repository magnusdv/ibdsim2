#' Segment distribution plot
#'
#' @param \dots One or several data frames, typically produced by
#'   [realisedInbreeding].
#' @param lab A character vector of labels used in the legend. Only relevant if
#'   `x` is an unnamed list of length > 1.
#' @param col A vector of the same length as `...`
#' @param alpha A transparency parameter for the scatter points.
#' @param ellipse A logical: Should confidence ellipses be added to the plot?
#' @param legendInside A logical controlling the legend placement.
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
#' # Simulations (increase N!))
#' s = lapply(peds, function(p)
#'   ibdsim(p, N = 10, ids = leaves(p), map = "uniform.sex.spec", verbose = FALSE))
#'
#' # Summarise autozygous regions
#' segs = lapply(s, realisedInbreeding)
#'
#' # Plot distributions
#' plotAutozygosity(segs, title = "Distribution of autozygous segments")
#'
#' @import ggplot2
#' @importFrom ribd inbreeding
#' @export
plotAutozygosity = function(..., lab = NULL, col = NULL, alpha = 1, 
                           ellipse = TRUE, title = NULL, 
                           xlab = "Segment count", 
                           ylab = "Average segment length (cM)",
                           legendInside = TRUE) {
  
  x = list(...)
  if(length(x) == 1 && !is.data.frame(x[[1]]))
    x = x[[1]]
  
  N = length(x)
  
  # Labels
  lab = lab %||% names(x) %||% as.character(1:N)
  if(N > 1 && any(lab == ""))
    stop2("Missing names in input list: ", lab)
  
  plotData = do.call(rbind, x)
  plotData$Relationship = factor(rep(lab, sapply(x, nrow)), levels = lab)
  
  # Theoretical expectation curves
  expected = sort(unique(plotData$fPed))
  genomeLen = plotData$genomeLen[1]
  
  # Remove 0 if included
  if(0 %in% expected) {
    message("No theoretical curve to draw for expectated value 0")
    expected = expected[expected > 0]
  }
  
  max.x = max(plotData$meanLen)
  max.y = max(plotData$nSeg)
  
  # Create plot
  g = ggplot(data = plotData, 
             aes_string(x = "meanLen", y = "nSeg", color = "Relationship")) + 
    geom_jitter(height = 0.35, alpha = alpha) +
    theme_bw(base_size = 15) + 
    scale_color_brewer(palette = "Set1") +
    labs(title = title, x = xlab, y = ylab) +
    {if(ellipse) stat_ellipse(size = 1.2)}
    
  # Add curves to plot
  if(length(expected)) {
    curveData = do.call(rbind, lapply(expected, function(v) {
      xmin = if(max.y > 0) genomeLen * v / max.y else 0
      xvec = seq(xmin, max.x, length = 50)
      data.frame(x = xvec, y = genomeLen * v / xvec, 
                 coeff = as.character(round(v, 3)))
    }))
  
  g = g + 
    geom_line(data = curveData, aes_string("x", "y", linetype = "coeff"), 
              lwd = 1, inherit.aes = FALSE) + 
    scale_linetype_manual(values = 1 + seq_along(expected), 
                          labels = paste0("f = ", expected),
                          name = "Expected")
  }
  
  # Fix legends
  g = g + 
    theme(legend.key.width = unit(0.9, "cm")) + 
    guides(color = if(N > 1) guide_legend(order = 1) else FALSE,
           linetype = guide_legend(order = 2, reverse = TRUE))
  
  if(legendInside) 
    g = g + theme(legend.position = c(.95, .95), 
                  legend.justification = c("right", "top"))
  
  g
}

 