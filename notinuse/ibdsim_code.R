### Code taken out from ibdsim()

... if(is.null(query)) return(simdata)

  coeffs = inbreeding(x); inbreds = which(coeffs > 0);
  if (length(inbreds) > 0) {
    inb = cbind(ID = inbreds, f = coeffs[inbreds]); rownames(inb) = rep("", length(inbreds))
    if (verbose) {
      cat("\nInbreeding coefficients:\n"); print(inb)
    }
  }

  runs <- lapply(simdata, function(h) sap.segments(h, sap = query))
  attr(runs, "total_map_length_Mb") = attr(simdata, "total_map_length_Mb")

  if (verbose) cat("\nResults:\n")
  stats = summary.ibd(runs, merged = merged, verbose = verbose)
  if (verbose) cat("\nTotal time used:", (proc.time() - starttime)[["elapsed"]], "seconds.\n")

  invisible(list(simdata = simdata, segments = runs, stats = stats))