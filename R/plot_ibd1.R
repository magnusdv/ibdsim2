#' Scatter plot of IBD=1 segment distributions
#'
#' Visualise and compare distributions of segments with IBD=1 between two
#' individuals.
#'
#' This function takes as input one or several complete outputs from the
#' [ibdsim()]. The underlying pedigree of each input is extracted, and
#' [pedtools::leaves()] is used to identify the two individuals to be
#' examined. (Hence it is required that each pedigree has exactly two leaves.)
#'
#' For each simulation the number of IBD=1 segments is plotted against the
#' average segment length. Finally, contour curves are added to plot,
#' corresponding to the theoretical (i.e., pedigree-based) value of
#' \eqn{\kappa_1} for each input relationship.
#'
#' @param ... One or several objects of class `genomeSimList`, i.e. outputs
#'   of [ibdsim()].
#' @param labels A character vector of the same length as above, containing
#'   descriptive plotting labels for the inputs.
#' @param alpha A transparency parameter for the scatter points.
#' @param ellipses A logical: Should confidence ellipses be added to the plot?
#' @param legend_inside A logical indicating wether the legend should be placed inside (default) or outside the plot window.
#' @param title,xlab,ylab Title and axis labels. Set to `NULL` to remove.
#' @return A ggplot2 plot object.
#'
#' @examples
#' ###
#' # Compare the distribution of IBD=1 segments
#' # between paternal and maternal half siblings.
#' ###
#' library(pedtools)
#' # Define the pedigrees
#' x.pat = halfSibPed()
#' x.mat = swapSex(x.pat, 1)
#'
#' # Simulate (increase 'sims'!)
#' map = "uniform.sex.spec"
#' sims = 10
#' s.pat = ibdsim(x.pat, sims = sims, map=map)
#' s.mat = ibdsim(x.mat, sims = sims, map=map)
#'
#' plot_ibd1(s.pat, s.mat, labels = c("HSpat", "HSmat"))
#'
#' @import ggplot2
#' @importFrom ribd ibd_kinship
#' @export
plot_ibd1 = function(..., labels, alpha = 1, ellipses = TRUE, legend_inside = FALSE, 
                     title = "Distribution of IBD segments", xlab = "Segment count", 
                     ylab = "Average segment length (cM)") {
  sims = list(...)
  if(missing(labels)) labels = seq_along(sims)
  stopifnot(length(labels) == length(sims))
  
  plot_data = lapply(seq_along(sims), function(i) {
    s = sims[[i]]
    ped = attr(s, 'pedigree')
    ids = leaves(ped)
    real = realised_kappa(s, id.pair=ids)
    if(any(real$Nsegments['Nseg2', ] > 0)) 
       message("Warning: Simulation list ", i, " includes IBD=2 segments. Expected kappa1 line will be wrong!")
    count = real$Nsegments["Nseg1", ]
    averlen = ifelse(count == 0, 0, real$kappa.realised['ibd1', ] * real$genomeLength / count)
    kinship_coeff = ibd_kinship(ped)[ids[1], ids[2]]
  
    data.frame(count = count,
               averlen = averlen, 
               genomeLength = real$genomeLength, 
               relation = labels[i], 
               ibd1_theory = kinship_coeff*4)
  })
  
  plot_data = do.call(rbind, plot_data)
  plot_data$relation = factor(plot_data$relation)
  
  max.x = max(plot_data$count)
  max.y = max(plot_data$averlen)
  
  g = ggplot(data=plot_data, aes_string(x="count", y="averlen", col="relation")) + 
    geom_jitter(width = 0.25, alpha = alpha) +
    theme_bw(base_size = 15) + 
    scale_color_manual(values = ggplotColors(nlevels(plot_data$relation))) +
    labs(title = title, x = xlab, y = ylab, col = "Relationship")
    
  if(ellipses) 
    g = g + stat_ellipse(aes_string(x = "count"), size=1.3)
  
  # Theoretical ibd1 curves
  ibd1_vals = sort(unique(plot_data$ibd1_theory))
  genomeLen = plot_data$genomeLength[1]
  
  curveData = do.call(rbind, lapply(ibd1_vals, function(ibd1) {
    xmin = genomeLen * ibd1 / max.y
    xvec = seq(xmin, max.x, length=10)
    data.frame(x = xvec, y = genomeLen*ibd1/xvec, kappa1 = as.character(round(ibd1, 3)))
  }))
  
  g = g + 
    geom_line(data = curveData, aes_string("x", "y", linetype = "kappa1"), lwd = 1, inherit.aes = F) + 
    scale_linetype_manual(values = 1 + seq_along(ibd1_vals), # avoid solid line
                          name = expression(Theoretical~kappa[1]))
  
  # Fix legends
  g = g + 
    theme(legend.key.width = unit(0.9, "cm")) + 
    guides(color = guide_legend(order = 1), 
           linetype = guide_legend(order = 2, reverse = T))
  
  if(legend_inside) 
    g = g + theme(legend.position = c(.95, .95), 
                  legend.justification = c("right", "top"))
  
  g
}

# Color generator
#' @importFrom grDevices hcl
ggplotColors = function(g){
  h = cumsum(c(15, rep(360/g, g-1)))
  g2 = ceiling(g/2)
  h = h[c(rbind(1:g2, 1:g2 + g2))[1:g]]
  hcl(h=h, c=100, l=50)
}
