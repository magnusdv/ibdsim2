library(ggplot2)
library(patchwork)

plotped = function(ped, ids, col, title, margin = 0.5) {
  fill = list(ids) |> setNames(col)
  suppressWarnings(
    plot(ped, fill = fill, title = title, cex.main = 1.5,
         cex = 1.25, autoscale = TRUE,  margin = margin)
  )
}

getSegmentData = function(sim, analysis, cutoff, unit) {
  if(is.null(sim))
    return(NULL)
  
  ids = extractIds(sim)
  
  # Extract segments according to analysis type
  pattern = setNames(list(ids), switch(analysis, Sharing = "carriers", Autozygosity = "autozygous"))
  segs = findPattern(sim, pattern, cutoff = cutoff, unit = unit)
  
  # Return stats
  segmentStats(segs, unit = unit, returnAll = TRUE)
}

generateIbdPlot = function(segData, analysis, cols, unit, observed = NULL) {
  # assume no entry or segData empty!
  labs = names(segData)
  
  perSimList = lapply(labs, function(lb) cbind(segData[[lb]]$perSim, Relationship = lb))
  perSim = do.call(rbind, perSimList)
  
  allSegsList = lapply(labs, function(lb) segData[[lb]]$allSegs)
  allSegs = data.frame(Length = unlist(allSegsList, use.names = FALSE), 
                       Relationship = rep(labs, lengths(allSegsList)))
  
  # Plot variables
  base = 14
  
  # Capitalise `unit`
  unit = switch(unit, cm = "cM", mb = "Mb", unit)
  
  # Plot1: Average length vs count
  g1 = ggplot(perSim, aes(.data$Count, .data$Average, color = .data$Relationship)) + 
    geom_jitter(width = 0.35, alpha = 0.5, shape = 1, show.legend = FALSE) +
    labs(x = "Number of segments", y = sprintf("Average segment (%s)", unit), col = NULL) +
    suppressWarnings(stat_ellipse(size = 1.2)) +
    theme_classic(base_size = base) + 
    theme(legend.position = "inside",
          legend.position.inside = c(.99, .99),
          legend.justification = c("right", "top"),
          legend.key.width = unit(1.4, "cm"),
          legend.text = element_text(size = 16), 
          plot.margin = margin(r = 40))
  
  # Plot 2: Total length distribution
  g2 = ggplot(perSim, aes(.data$Total, color = .data$Relationship)) + 
    geom_density(aes(fill = .data$Relationship), alpha = 0.1, size = 1.2, show.legend = FALSE) +
    labs(x = switch(analysis, 
                    Sharing = sprintf("Total length IBD (%s)", unit), 
                    Autozygosity = sprintf("Total autozygosity (%s)", unit)),
         y = "density") + 
    scale_y_continuous(expand = expansion(mult = c(0,0.04))) + 
    theme_classic(base_size = base) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_text(size = 12))
  
  
  # Plot 3: Count distribution
  g3 = ggplot(perSim, aes(.data$Count, color = .data$Relationship)) + 
    geom_density(aes(fill = .data$Relationship), alpha = 0.1, size = 1.2, show.legend = FALSE) +
    labs(x = "Number of segments", y = "density") +
    scale_y_continuous(expand = expansion(mult = c(0,0.04))) + 
    theme_classic(base_size = base) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_text(size = 12))
  
  # Plot 4: Segment length distribution
  g4 = ggplot(allSegs, aes(.data$Length, color = .data$Relationship)) + 
  geom_density(aes(fill = .data$Relationship), alpha = 0.1, size = 1.2, show.legend = FALSE) +
  labs(x = sprintf("Individual segment length (%s)", unit), y = "density") +
  coord_cartesian(xlim = c(0, quantile(allSegs$Length, 0.999))) + # avoid very long tails!
  
  scale_y_continuous(expand = expansion(mult = c(0,0.04))) + 
  theme_classic(base_size = base) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 12))

  # Add observed data
  if(!is.null(observed)) {
    g1 = g1 + annotate("point", x = observed$nseg, observed$mean, size = 3.5, stroke = 2, shape = 4, colour = 1)
    g2 = g2 + annotate("text", x = observed$total, y = 0, label = "\u25B2", vjust = 0, size = 6)
    g3 = g3 + annotate("text", x = observed$nseg, y = 0, label = "\u25B2", vjust = 0, size = 6)
    if(length(observed$lengths))
      g4 = g4 + annotate("text", x = observed$lengths, y = 0, label = "|", vjust = 0, size = 3)
  }
  
  (g1 | (g2 / g3 / g4)) & # plot_layout(guides = 'collect') & 
  scale_color_manual(values = cols)
}


calculateMargin = function(w) {
  side = if(w < 200) 0.5 else if(w < 300) 1 else if(w < 400) 2 else 3
  c(0.5, side, 3, side)
}
