#' Scatter plot of IBD=1 segment distributions
#'
#' Visualise and compare distributions of segments with IBD=1 between two
#' individuals.
#'
#' This function takes as input one or several complete outputs from the
#' [ibdsim()], and investigates the IBD sharing of two indicated individuals
#' from each. (See the `ids` entry in the parameter list for how the individuals
#' are indicated.)
#'
#' For each simulation the number of IBD=1 segments is plotted against the
#' average segment length. Finally, contour curves are added to plot,
#' corresponding to the theoretical (i.e., pedigree-based) value of
#' \eqn{\kappa_1} for each input relationship.
#'
#' @param ... One or several objects of class `genomeSimList`, i.e. outputs of
#'   [ibdsim()].
#' @param labels A character vector of the same length as `...`, containing
#'   descriptive plotting labels for the inputs.
#' @param pairs A list of the same length as `...`, containing a pair of ID labels for
#'   each input pedigree. The list names must equal (as a set) the `labels`
#'   argument. Two shortcuts are possible: If `pairs` is a vector of two ID
#'   labels, these will be used for all pedigrees. Finally, if `pairs` is the word
#'   "leaves", then `pedtools::leaves()` is used to extract the ID labels of
#'   each pedigree (they must all have exactly two leaves for this to work).
#' @param alpha A transparency parameter for the scatter points.
#' @param ellipses A logical: Should confidence ellipses be added to the plot?
#' @param legend_inside A logical indicating wether the legend should be placed
#'   inside (default) or outside the plot window.
#' @param title,xlab,ylab Title and axis labels. Set to `NULL` to remove.
#' @return A ggplot2 plot object.
#'
#' @examples
#' ###
#' # Example 1: Comparison of IBD segment distributions
#' # between paternal and maternal half siblings.
#' ###
#'
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
#' # By default, the IBD segments of the "leaves" are computed and plotted
#' plot_ibd1(s.pat, s.mat, labels = c("HSpat", "HSmat"))
#'
#' ###
#' # Example 2: Half siblings vs grandparent/grandchild
#' ###
#'
#' # Only one pedigree needed here
#' y = addSon(x.pat, 5)
#' 
#' # Simulate 
#' s = ibdsim(y, sims = sims, map=map)
#'
#' # Indicate the pairs explicitly this time.
#' # List names are used as labels in the plot
#' plot_ibd1(s, pairs = list(HS = 4:5, GR = c(1,7)))
#' 
#' 
#' @import ggplot2
#' @importFrom ribd kinship
#' @export
plot_ibd1 = function(..., labels, pairs = "leaves", alpha = 1, ellipses = TRUE, 
                     legend_inside = TRUE, title = "Distribution of IBD segments", 
                     xlab = "Segment count", ylab = "Average segment length (cM)") {
  sims = list(...)
  N = length(sims)
  if(N == 0) 
    stop2("The first argument must be a `genomeSimList` object, typically an output of ibdsim()")
  
  # If a single sim object and multiple pairs: Repeat sim 
  if(N == 1 && is.list(pairs))
    sims = rep(sims, N <- length(pairs))
  
  # Check labels vs names(pairs)
  if(!missing(labels) && !is.null(names(pairs)))
    stopifnot(length(labels) == N, !anyDuplicated.default(labels), 
              length(names(pairs)) == N, setequal(labels, names(pairs)))
  
  # Missing "labels" argument?
  if(missing(labels) || is.null(labels))
      labels = names(pairs) %||% 1:N
  else {
    # If `pairs` has names, use them to sort according to labels
    if(!is.null(names(pairs)))
      pairs = pairs[labels]
  }
  
  # Shortcut options for the `pairs` parameter
  if(identical(pairs, "leaves")) { 
    pairs = lapply(sims, function(s) leaves(attr(s, 'pedigree')))
  }
  else if(is.atomic(pairs)) { # Same pair for all peds
    pairs = rep(list(pairs), N)
  }
  
  # Check that the individuals in `pairs` exist
  for(i in 1:N) {
    ped = attr(sims[[i]], 'pedigree')
    ids = pairs[[i]]
    lab = labels[i]
    if(length(ids) != 2 || ids[1] == ids[2]) 
      stop2("The `pairs` entry for pedigree ", lab, " is not a valid pair: ", ids)
    int_ids = match(ids, labels(ped))
    if (anyNA(int_ids))
      stop2("Unknown ID label in pedigree ", lab, ": ", ids[is.na(int_ids)])
  }
  
  plot_data = lapply(1:N, function(i) {
    s = sims[[i]]
    ped = attr(s, 'pedigree')
    ids = pairs[[i]]
    real = realised_kappa(s, id.pair=ids)
    if(any(real$Nsegments['Nseg2', ] > 0)) 
       message("Warning: Simulation list ", i, " includes IBD = 2 segments. Expected 'kappa_1 curve' will be wrong!")
    count = real$Nsegments["Nseg1", ]
    averlen = ifelse(count == 0, 0, real$kappa.realised['ibd1', ] * real$genomeLength / count)
    kinship_coeff = kinship(ped)[ids[1], ids[2]]
  
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
  
  # Create plot
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
  
  # Remove 0 if included
  if(0 %in% ibd1_vals) {
    message("No theoretical curve to draw for kappa_1 = 0")
    ibd1_vals = ibd1_vals[ibd1_vals > 0]
  }
  
  # Add curves to plot
  if(length(ibd1_vals) > 0) {
    curveData = do.call(rbind, lapply(ibd1_vals, function(ibd1) {
      xmin = if(max.y > 0) genomeLen * ibd1 / max.y else 0
      xvec = seq(xmin, max.x, length = 10)
      data.frame(x = xvec, y = genomeLen*ibd1/xvec, kappa1 = as.character(round(ibd1, 3)))
    }))
    
    g = g + 
      geom_line(data = curveData, aes_string("x", "y", linetype = "kappa1"), 
                lwd = 1, inherit.aes = F) + 
      scale_linetype_manual(values = 1 + seq_along(ibd1_vals), # avoid solid line
                            name = expression(Theoretical~kappa[1]))
  }
  
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
  hcl(h = h, c = 100, l = 50)
}
