#' IBD simulation
#'
#' This is the main function of the package. Gene dropping and recombination of
#' chromosomes is simulated down the pedigree, either unconditionally or
#' conditional on given allele patterns.
#'
#' Each simulation starts by unique alleles being distributed to the pedigree
#' founders. In each subsequent meiosis, homologue chromosomes are made to
#' recombine according to a renewal process along the four-strand bundle, with
#' chi square distributed waiting times. (For comparison purposes, Haldane's
#' Poisson model for recombination is also implemented.)
#'
#' Recombination rates are sex-dependent, and vary along each chromosome
#' according to the recombination map specified by the \code{map} parameter. By
#' default, the complete Decode map of the human autosome is used (see
#' References). If \code{map="uniform.sex.spec"}, the genetic chromosome
#' \emph{lengths} are as in the Decode map, but the recombination rate is kept
#' constant along each chromosome. If \code{map="uniform.sex.aver"}, sex
#' averaged genetic chromosome lengths are used (and constant recombination
#' rates along each chromosome).
#'
#' IBD patterns are described as combinations of Single Allele Patterns (SAPs).
#' A SAP is a specification for a given allele of the number of copies carried
#' by various individuals, and must be given as a list of numerical vectors
#' (containing ID labels) named '0', '1', '2', 'atleast1' and 'atmost1' (some of
#' these can be absent or NULL; see Examples).
#'
#' If a condition SAP is given (i.e. if \code{condition} is non-null),
#' simulation of each complete chromosome set (all autosomes by default) is
#' performed as follows: A 'disease chromosome' is sampled at random (using the
#' physical chromosome lengths as weights), followed by a random 'disease locus'
#' on this chromosome.  For this chromosome, gene dropping down the pedigree is
#' carried out in such a way that the 'disease locus' has the condition SAP. (In
#' a bit more detail: First, the program computes all possible sets of obligate
#' carriers, with suitable weights, and samples one of these.  Included in the
#' obligate carriers will e exactly one founder, one of whose alleles is taken
#' as the 'disease allele'. In each meiosis involving obligate carriers,
#' recombination is performed as usual, but the strand carrying the 'disease
#' allele' is always the one passed on.) For the other chromosomes, simulation
#' is done unconditionally.
#'
#' @param x A pedigree in the form of a \code{\link[pedtools]{ped}} object.
#' @param sims A positive integer indicating the number of simulations.
#' @param condition A single allele pattern (SAP), i.e., a list with
#'   numerical entries named "0", "1", "2", "atleast1", "atmost1".
#' @param map The genetic map(s) to be used in the simulations: One of the
#'   character strings "decode", "uniform.sex.spec", "uniform.sex.aver". (See
#'   Details.)
#' @param chromosomes A numeric vector indicating chromosome numbers, or either
#'   of the words "AUTOSOMAL" or "X". The default is 1:22, i.e. the human
#'   autosomes.
#' @param model A character indicating the statistical model for recombination:
#'   Either "chi" (the default) or "haldane". (See details.)
#' @param skip.recomb A numeric containing individuals whose meioses should be
#'   simulated without recombination (i.e. a random strand is passed on to each
#'   offspring). If NULL, nobody is skipped. The default value (the character
#'   "noninf_founders") computes the set of pedigree founders that cannot be
#'   carriers of the alleles described in the \code{condition} SAPs.
#' @param seed An integer to be passed on to \code{\link{set.seed}}).
#' @param verbose A logical.
#'
#' @return The simulated genomes are invisibly returned.
#'
#' @examples
#'
#' z = pedtools::halfCousinsPed(0)
#' plot(z)
#' res = ibdsim(z, sims=1)
#' res
#'
#' @export
#' @importFrom pedtools is.ped has_parents_before_children parents_before_children
#' @importFrom assertthat assert_that is.count  
ibdsim = function(x, sims, condition=NULL, map="decode", chromosomes=NULL,
                  model="chi", skip.recomb = "noninf_founders", seed=NULL, verbose=TRUE) {
  # Check input
  assert_that(is.ped(x), is.count(sims))
  
  if(!model %in% c("chi", "haldane"))
    stop("Argument 'model' must be either 'chi' or 'haldane'.")
  model_string = if(model=="chi") "Chi square renewal process" else "Haldane's poisson process"    
  
  # Ensure that parents preceede their children
  if (!pedtools::has_parents_before_children(x)) {
    message("Reordering so that all parents preceede their children")
    x = pedtools::parents_before_children(x)
  }
  
  # Start timer
  starttime = proc.time()

  # Load map and extract chromosome names.
  map = loadMap(map, chrom = chromosomes)
  mapchrom = attr(map, "chromosome") # is NULL if map contains several
  if (is.null(mapchrom)) 
    mapchrom = sapply(map, attr, "chromosome")

  if (verbose) {
    cond_str = if (is.null(condition)) "unconditional" else "conditional"
    chroms = paste(mapchrom, collapse=",")
    
    print(glue::glue("
    Performing {cond_str} simulation.
    Chromosomes: {chroms}
    Recombination model: {model_string}
    Number of simulations: {sims}"))
  }
  
  # Seed for random sampling
  set.seed(seed)
  
  # Setup for conditional simulation
  dischr = numeric(sims)
  if (!is.null(condition)) {
    
    # Sample chromosome carrying the conditional locus
    if (length(map) == 1) 
      dischr[] = rep.int(mapchrom, sims) 
    else {
      chromlengths = sapply(map, attr, "length_Mb")
      dischr[] = sample(mapchrom, size = sims, replace = T, prob = chromlengths)
    }
    # Process condition SAP
    oblig.saps = sample.obligates(x, condition, sims)
  }
  
  # Determine ped members where recombination should be skipped.
  if (!is.null(skip.recomb)) {
    if (skip.recomb == "noninf_founders") {
      cafs = x$FOUNDERS
      if (!is.null(condition)) 
        cafs = intersect(cafs, .CAFs(x, condition))
      skip.recomb = setdiff(x$FOUNDERS, cafs)
    }
    if (length(skip.recomb) > 0 && verbose) 
      message("Skipping recombination in:", paste(skip.recomb, collapse = ","))
  }
  
  # The actual simulations: One sim at the time; each chromosome in turn 
  genomeSimList = lapply(1:sims, function(i) {
    lapply(map, function(m) {
      # TODO: should dischr[sims] below be dischr[s]?
      cond = if (dischr[sims] == attr(m, "chromosome")) oblig.saps[[i]] else NULL 
      genedrop(x, map = m, condition = cond, model = model, skip.recomb = skip.recomb)
    })
  })
  
  if (verbose) {
    elapsed = (proc.time() - starttime)[["elapsed"]]
    message("Simulation finished in ", elapsed, " seconds.")
  }
  
  # Various attributes of the simulation call
  attribs = list(pedigree = x, 
                 skipped = skip.recomb,
                 condition = condition,
                 genome_length_Mb = attr(map, "length_Mb"),
                 chromosomes = mapchrom,
                 model = model_string)
  
  # Add attributes and class to each genomeSim
  genomeSimList = lapply(genomeSimList, `attributes<-`, c(attribs, class="genomeSim"))
  
  # Add attributes and class to the entire list 
  attributes(genomeSimList) = c(attribs, class = "genomeSimList")
  
  genomeSimList
}


sample.obligates = function(x, condition, sims) {
  obligate_ones = obligate.carriers(x, condition)
  complete.saps = lapply(obligate_ones, function(ones) {
    sap = condition; sap[["1"]] = ones; sap
  })
  if (length(complete.saps) == 1) {
    cat("For the disease chromosome I'm conditioning on the following SAP:\n"); .printSAP(complete.saps[[1]])
  }
  else {
    cat("For the disease chromosome I'm sampling condition SAPs among the following:\n")
    for (i in 1:length(complete.saps)) {
      cat("SAP ", i, ":\n", sep = ""); .printSAP(complete.saps[[i]])
    }
  }
  weight = sapply(obligate_ones, function(vec) .5^(length(vec) - 1))
  oblig.samples = sample(complete.saps, size = sims, replace = TRUE, prob = weight)
}

