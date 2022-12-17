#' Scatter plots of IBD segment distributions
#'
#' Visualise and compare count/length distributions of IBD segments. Two types
#' are currently implemented: Segments of autozygosity (for a single person) and
#' segments with (pairwise) IBD state 1.
#'
#' This function takes as input one or several complete outputs from the
#' [ibdsim()], and produces a scatter plot of the number and average length of
#' IBD segments from each.
#'
#' Contour curves are added to plot, corresponding to the
#' theoretical/pedigree-based values: either inbreeding coefficients (if `type =
#' "autozygosity"`) or \eqn{\kappa_1} (if `type = "ibd1"`).
#'
#' @param \dots One or several objects of class `genomeSimList`, typically
#'   created by [ibdsim()]. They can be entered separately or as a `list`.
#' @param type A string indicating which segments should be plotted. Currently,
#'   the allowed entries are "autozygosity" and "ibd1".
#' @param ids A list of the same length as `...`, where each entry contains one
#'   or two ID labels (depending on `type`). By default (NULL), these labels are
#'   extracted from the inputs in `...`.
#'
#'   Two other short-cuts are possible: If a single vector is given, it is
#'   repeated for all pedigrees. Finally, if `ids` is the word "leaves" then
#'   `pedtools::leaves()` is used to extract labels in each pedigree.
#' @param labels An optional character vector of labels used in the legend. If
#'   NULL, the labels are taken from `names(...)`.
#' @param unit Length unit; either "cm" (centiMorgan) or "mb" (megabases).
#' @param col An optional colour vector of the same length as `...`.
#' @param shape A vector with point shapes, of the same length as `...`.
#' @param alpha A transparency parameter for the scatter points.
#' @param ellipses A logical: Should confidence ellipses be added to the plot?
#' @param legendInside A logical controlling the legend placement.
#' @param title,xlab,ylab Title and axis labels.
#'
#' @examples
#'
#' # Simulation parameters used in the below examples.
#' map = uniformMap(M = 10)   # recombination map
#' N = 5                      # number of sims
#'
#' # For more realistic results, replace with e.g.:
#' # map = loadMap("decode19")
#' # N = 1000
#'
#'
#' #################################################################
#' # EXAMPLE 1
#' # Comparison of IBD segment distributions
#' # between paternal and maternal half siblings.
#' #################################################################
#'
#' # Define the pedigrees
#' xPat = halfSibPed()
#' xMat = swapSex(xPat, 1)
#'
#' simPat = ibdsim(xPat, N = N, map = map)
#' simMat = ibdsim(xMat, N = N, map = map)
#'
#' # By default, the IBD segments of the "leaves" are computed and plotted
#' plotSegmentDistribution(simPat, simMat, type = "ibd1", ids = 4:5,
#'                         labels = c("HSpat", "HSmat"))
#'
#' #################################################################
#' # EXAMPLE 2
#' # Half siblings vs half uncle vs grandparent/grandchild
#' #################################################################
#'
#' # Only one pedigree needed here
#' x = addSon(halfSibPed(), 5)
#'
#' s = ibdsim(x, N = N, map = map)
#'
#' # Indicate the pairs explicitly this time.
#' ids = list(GR = c(2,7), HS = 4:5, HU = c(4,7))
#'
#' # List names are used as labels in the plot
#' plotSegmentDistribution(s, type = "ibd1", ids = ids, shape = 1:3)
#'
#'
#' #################################################################
#' # EXAMPLE 3
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
#' # Simulations
#' s = lapply(peds, function(p)
#'   ibdsim(p, N = N, ids = leaves(p), verbose = FALSE, map = map))
#'
#' # Plot distributions
#' plotSegmentDistribution(s, type = "autoz", title = "Autozygous segments")
#'
#' @export
plotSegmentDistribution = function(..., type = c("autozygosity", "ibd1"), 
                                   ids = NULL, unit = "cm", labels = NULL, col = NULL, 
                                   shape = 1, alpha = 1, ellipses = TRUE, 
                                   title = NULL, xlab = NULL, ylab = NULL,
                                   legendInside = TRUE) {

  sims = list(...)
  N = length(sims)
  if(N == 0 || !(inherits(sims[[1]], "genomeSimList") || inherits(sims[[1]][[1]], "genomeSimList")))
    stop2("The first argument must be a `genomeSimList` object or a list of such")
  
  # If input was already a list, flatten it
  if(N == 1 && !inherits(sims[[1]], "genomeSimList")) {
    sims = sims[[1]]
    N = length(sims)
  }
  
  # Check that all ... are ok
  if(!all(ok <- sapply(sims, inherits, "genomeSimList")))
    stop2("Bad input detected: ", setdiff(names(sims)[!ok], ""), 
          "\n(This function does not allow abbreviated argument names)")
  
  # If a single sim object and multiple pairs: Repeat sim 
  if(N == 1 && is.list(ids) && length(ids) > 1) {
    N = length(ids)
    sims = rep(sims, N)
  }
  
  # If both sims and ids have names, use them to sort ids
  simnms = names(sims)
  idsnms = names(ids)
  if(is.list(ids) && !is.null(idsnms) && !is.null(simnms)) {
    stopifnot(setequal(simnms, idsnms))
    ids = ids[simnms]
  }
  
  # Other options for the `ids` parameter
  if(is.null(ids))
    ids = lapply(sims, extractIds)
  else if(identical(ids, "leaves"))
    ids = lapply(sims, function(s) leaves(attr(s, 'pedigree')))
  else if(is.atomic(ids))
    ids = rep(list(ids), N)
  
  # Check `ids` compatibility
  for(i in 1:N) {
    ped = attr(sims[[i]], 'pedigree')
    if(!all(ids[[i]] %in% labels(ped)))
      stop2("Unknown ID label in pedigree ", i, ": ", setdiff(ids[[i]], labels(ped)))
  }
  
  # Set labels as list names
  labs = labels %||% simnms %||% names(ids) %||% as.character(1:N)
  stopifnot(length(labs) == N, !anyDuplicated.default(labs))
  names(sims) = names(ids) = labs
  
  col = if(is.null(col)) (1 + 1:N) else rep(col, length.out = N)
  shape = rep(shape, N)
  xlab = xlab %||% "Number of segments"
  ylab = ylab %||% sprintf("Average length (%s)", unit) 
  
  PLOTFUN = switch(match.arg(type), 
                   autozygosity = plotSegmentDistribution.autoz, 
                   ibd1 = plotSegmentDistribution.ibd1)
  
  PLOTFUN(sims, ids, unit = unit, col = col, shape = shape, alpha = alpha, ellipses = ellipses, 
          title = title, xlab = xlab, ylab = ylab, legendInside = legendInside)
}
  

#' @importFrom ribd inbreeding 
plotSegmentDistribution.autoz = function(sims, ids, unit, col = NULL, shape = 1, alpha = 1, ellipses = TRUE, 
                                         title = NULL, xlab = NULL, ylab = NULL, legendInside = TRUE) {
  
  N = length(sims)
  if(any(lengths(ids) != 1))
    for(i in 1:N) if(length(ids[[i]]) != 1) 
      stop2("The `ids` entry for pedigree ", i, " is not a single individual: ", ids[[i]])
  
  labs = names(sims)
  
  plotDatList = lapply(1:N, function(i) {
    s = sims[[i]]
    id = ids[[i]]
    real = realisedInbreeding(s, id = id, unit = unit)$perSimulation
    cbind(real[c('nSeg', 'meanLen')], relation = labs[i])
  })
  
  plotDat = do.call(rbind, plotDatList)
  plotDat$relation = factor(plotDat$relation, levels = unique.default(plotDat$relation))
  
  ### Theoretical expectation curves
  fPed = sapply(1:N, function(i) inbreeding(attr(sims[[i]], "pedigree"), ids = ids[[i]]))
  
  expect.args = list(values = fPed, 
                     physRange = .getLen(sims[[1]], unit),
                     label = expression(Expected~italic(f)))
  
  # Create the plot
  .plotSegDist(plotDat, col = col, shape = shape, alpha = alpha, ellipses = ellipses, title = title, xlab = xlab, 
               ylab = ylab, expect.args = expect.args, legendInside = legendInside)
  
}


#' @importFrom ribd kinship
plotSegmentDistribution.ibd1 = function(sims, ids, unit, col = NULL, shape = 1, alpha = 1, ellipses = TRUE, 
                                        title = NULL, xlab = NULL, ylab = NULL, legendInside = TRUE) {
  
  N = length(sims)
  if(any(lengths(ids) != 2)) {
    for(i in 1:N) if(length(ids[[i]]) != 2) 
      stop2("The `ids` entry for pedigree ", i, " is not a valid pair: ", ids[[i]])
  }
  
  labs = names(sims)
  
  plotDatList = lapply(1:N, function(i) {
    s = sims[[i]]
    ids = ids[[i]]
    real = realisedKappa(s, ids = ids, unit = unit)$perSimulation
    
    if(any(real$nSeg2 > 0)) 
      message("Warning: Simulation list ", i, " includes IBD = 2 segments. Expected 'kappa_1 curve' will be wrong!")
    
    L = .getLen(s, unit)
    nSeg = real$nSeg1
    meanLen = ifelse(nSeg > 0, real$k1 * L / nSeg, 0)
    data.frame(nSeg = nSeg,  meanLen = meanLen, relation = labs[i])
  })
  
  plotDat = do.call(rbind, plotDatList)
  plotDat$relation = factor(plotDat$relation, levels = unique.default(plotDat$relation))
  
  ### Theoretical expectation curves
  kappa1 = sapply(1:N, function(i) 
    4 * kinship(attr(sims[[i]], "pedigree"), as.character(ids[[i]])))
  
  expect.args = list(values = kappa1, 
                     physRange = .getLen(sims[[1]], unit),
                     label = expression(Expected~kappa[1]))
  
  # Create the plot
  .plotSegDist(plotDat, col = col, shape = shape, alpha = alpha, ellipses = ellipses, title = title, xlab = xlab, 
               ylab = ylab, expect.args = expect.args, legendInside = legendInside)
  
}

#' @import ggplot2
.plotSegDist = function(plotDat, col, shape, alpha, ellipses, title, xlab, ylab, expect.args, legendInside = TRUE) {
  nRel = nlevels(plotDat$relation)
  
  g = ggplot(plotDat, aes_string(x = "nSeg", y = "meanLen", color = "relation", shape = "relation")) + 
    geom_jitter(width = 0.35, alpha = alpha) +
    labs(title = title, x = xlab, y = ylab, col = "Relationship", shape = "Relationship") +
    theme_bw(base_size = 15) + 
    scale_color_manual(values = col) +
    scale_shape_manual(values = rep(shape, length.out = nRel)) +
    guides(color = if(nRel > 1) guide_legend(order = 1, override.aes = list(size = 2, alpha = 1)) 
                   else "none",
           shape = if(nRel > 1) guide_legend(order = 1) 
                   else "none") + 
    theme(legend.key.width = unit(0.9, "cm"))
  
  if(legendInside)
    g = g + theme(legend.position = c(.95, .95),
                  legend.justification = c("right", "top"))
  
  if(ellipses) 
    g = g + stat_ellipse(size = 1.2, show.legend = FALSE)
  
  #### Expectation curves ###
  vals = expect.args$values
  vals = unique.default(setdiff(vals, 0))
  
  # Data for expectation curves
  if(length(vals)) {
    max.x = max(plotDat$nSeg)
    max.y = max(plotDat$meanLen)
    L = expect.args$physRange
    
    curveDat = do.call(rbind, lapply(vals, function(v) {
      xmin = if(max.y > 0) L * v / max.y else 0
      xvec = seq(xmin, max.x, length.out = 50)
      data.frame(x = xvec, y = L * v / xvec, 
                 coeff = as.character(round(v, 3)), stringsAsFactors = TRUE)
    }))
    
    g = g + 
      geom_line(data = curveDat, 
                aes_string("x", "y", linetype = "coeff"), 
                lwd = 1, inherit.aes = FALSE) + 
      scale_linetype_manual(values = 1 + seq_along(vals), 
                            name = expect.args$label)
  }
  
  g
}

.getLen = function(s, unit) {
  if(unit == "mb")
    return(attr(s, "physRange"))
  
  mapl = attr(s, "mapLen")
  if(isXsim(s))
    mapl[["female"]]
  else 
    mean(mapl)
}
