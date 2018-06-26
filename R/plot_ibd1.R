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
#'
#' @return A ggplot2 plot object.
#'
#' @examples
#' ###
#' # Compare the distribution of IBD=1 segments
#' # between paternal and maternal half siblings.
#' ###
#'
#' # Define the pedigrees
#' x.pat = pedtools::halfCousinsPed(0)
#' x.mat = pedtools::swapSex(x.pat, 1)
#'
#' # Simulate (increase 'sims'!)
#' map = "uniform.sex.spec"
#' s.pat = ibdsim(x.pat, sims = 10, map=map)
#' s.mat = ibdsim(x.mat, sims = 10, map=map)
#'
#' plot_ibd1(s.pat, s.mat, labels=c("HSpat", "HSmat"))
#'
#' @import ggplot2
#' @importFrom kinship2 kinship
#' @importFrom pedtools leaves as.kinship2_pedigree
#' @export
plot_ibd1 = function(..., labels, alpha=1) {
  sims = list(...)
  if(missing(labels)) labels = seq_along(sims)
  stopifnot(length(labels) == length(sims))
  
  plot_data = lapply(seq_along(sims), function(i) {
    s = sims[[i]]
    ped = attr(s, 'pedigree')
    ids = pedtools::leaves(ped)
    real = realised_kappa(s, id.pair=ids)
    if(any(real$Nsegments['Nseg2', ] > 0)) 
       message("Warning: Simulation list ", i, " includes IBD=2 segments. Expected kappa1 line will be wrong!")
    count = real$Nsegments["Nseg1", ]
    averlen = ifelse(count == 0, 0, real$kappa.realised['ibd1', ] * real$genomeLength / count)
    kinship_coeff = kinship2::kinship(as.kinship2_pedigree(ped))[ids[1], ids[2]]
    
    data.frame(count = count, 
               averlen = averlen, 
               genomeLength = real$genomeLength, 
               relation = labels[i], 
               ibd1_theory = kinship_coeff*4)
  })
  plot_data = do.call(rbind, plot_data)
  plot_data$relation = factor(plot_data$relation)
  
  range.x = range(plot_data$count)
  max.x = max(plot_data$count)
  max.y = max(plot_data$averlen)
  
  g = ggplot(data=plot_data, aes_string("count", "averlen", col="relation")) + 
    theme_bw(base_size=14) + 
    scale_color_manual(values = ggplotColors(nlevels(plot_data$relation))) +
    geom_jitter(width=0.4, height=0.5, alpha=alpha) + 
    stat_ellipse(size=1.3) + 
    labs(x="Segment count", y = "Average length (cM)", title="Distribution of IBD segments")
  
  # Theoretical ibd1 curves
  ibd1_vals = unique(plot_data$ibd1_theory)
  genomeLen = plot_data$genomeLength[1]
  
  pointlist = lapply(ibd1_vals, function(ibd1) {
    xmin = genomeLen * ibd1 / max.y
    xvec = seq(xmin, max.x, length=100)
    data.frame(IBD1=ibd1, x=xvec, y=genomeLen*ibd1/xvec)
  })
  
  contour.data = do.call(rbind, pointlist)
  contour.data$IBD1 = as.factor(contour.data$IBD1)
  
  ymin = sapply(pointlist, function(df) df$y[nrow(df)])
  kappa_labels = paste("kappa[1] == ", ibd1_vals)
  
  g = g + 
    geom_line(data=contour.data, aes_string("x", "y", group="IBD1"), lty=2, lwd=1, inherit.aes = F) +
    annotate("text", max.x + 1, ymin, label=kappa_labels, hjust="middle", vjust="top", parse=T)
  
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
