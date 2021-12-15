#' IBD simulation
#'
#' This is the main function of the package, simulating the recombination
#' process in each meioses of a pedigree. The output summarises the IBD segments
#' between all or a subset of individuals.
#'
#' Each simulation starts by unique alleles (labelled 1, 2, ...) being
#' distributed to the pedigree founders. In each meiosis, homologue chromosomes
#' are made to recombine according to the value of `model`:
#'
#' * `model = "haldane"`: In this model, crossover events are modelled as a
#' Poisson process along each chromosome.
#'
#' * `model = "chi"` (default): This uses a renewal process along the
#' four-strand bundle, with waiting times following a chi square distribution.
#'
#' Recombination rates along each chromosome are determined by the `map`
#' parameter. The default value ("decode19") loads a thinned version of the
#' recombination map of the human genome published by Halldorsson et al (2019).
#'
#' In many applications, the fine-scale default map is not necessary, and should
#' be replaced by simpler maps with constant recombination rates. See
#' [uniformMap()] and [loadMap()] for ways to produce such maps.
#'
#' @param x A [pedtools::ped()] object.
#' @param N A positive integer indicating the number of simulations.
#' @param ids A subset of pedigree members whose IBD sharing should be analysed.
#'   If NULL, all members are included.
#' @param map The genetic map to be used in the simulations: Allowed values are:
#'
#'   * a `genomeMap` object, typically produced by [loadMap()]
#'
#'   * a single `chromMap` object, for instance as produced by [uniformMap()]
#'
#'   * a character, which is passed on to [loadMap()] with default parameters.
#'   Currently the only valid option is "decode19" (or abbreviations of this).
#'
#'   Default: "decode19".
#'
#' @param model Either "chi" or "haldane", indicating the statistical model for
#'   recombination (see details). Default: "chi".
#' @param skipRecomb A vector of ID labels indicating individuals whose meioses
#'   should be simulated without recombination. (Each child will then receive a
#'   random strand of each chromosome.) The default action is to skip
#'   recombination in founders who are uninformative for IBD sharing in the
#'   `ids` individuals.
#'
#' @param seed An integer to be passed on to [set.seed()]).
#' @param verbose A logical.
#'
#' @return A list of `N` objects of class `genomeSim`.
#'
#'   A `genomeSim` object is essentially a numerical matrix describing the
#'   allele flow through the pedigree in a single simulated. Each row
#'   corresponds to a chromosomal segment. The first 4 columns describe the
#'   segment (chromosome, start, end, length), and are followed by two columns
#'   (paternal allele, maternal allele) for each of the `ids` individuals.
#'
#'   If `ids` has length 1, a column named "Aut" is added, whose entries are 1
#'   for autozygous segments and 0 otherwise.
#'
#'   If `ids` has length 2, two columns are added:
#'
#'   * `IBD` : The IBD status of each segment (= number of alleles shared
#'   identical by descent). For a given segment, the IBD status is either 0, 1,
#'   2 or NA. If either individual is inbred, they may be autozygous in a
#'   segment, in which case the IBD status is reported as NA. With inbred
#'   individuals the `Sigma` column (see below) is more informative than the
#'   `IBD` column.
#'
#'   * `Sigma` : The condensed identity ("Jacquard") state of each segment,
#'   given as an integer in the range 1-9. The numbers correspond to the
#'   standard ordering of the condensed states. In particular, for non-inbred
#'   individuals the states 9, 8, 7 correspond to IBD status 0, 1, 2
#'   respectively.
#'
#' @references Halldorsson et al. _Characterizing mutagenic effects of
#'   recombination through a sequence-level genetic map._ Science 363, no. 6425
#'   (2019).
#'
#' @examples
#'
#' hs = halfSibPed()
#' ibdsim(hs, N = 2, map = uniformMap(M = 1), ids = 4:5)
#'
#' # Full sib mating: all 9 states are possible
#' x = fullSibMating(1)
#' sim = ibdsim(x, N = 1, ids = 5:6, map = uniformMap(M = 10), seed = 1)
#' s = sim[[1]]
#' stopifnot(setequal(s[, 'Sigma'], 1:9))
#'
#' @export
ibdsim = function(x, N = 1, ids = labels(x), map = "decode",
                  model = c("chi", "haldane"), skipRecomb = NULL, 
                  seed = NULL, verbose = TRUE) {
  # Check input
  if(!is.ped(x))
    stop2("The first argument must be a `ped` object")
  if(!isCount(N))
    stop2("`N` must be a positive integer")
  if(!all(founderInbreeding(x) %in% c(0,1)))
    stop2("Founder inbreeding coefficients other than 0 and 1 are not allowed")
  
  # Ensure that parents precede their children
  if (!hasParentsBeforeChildren(x)) {
    if(verbose) message("Reordering so that all parents precede their children")
    x = parentsBeforeChildren(x)
  }
  
  # Start timer
  starttime = Sys.time()

  # Load map and extract chromosome names
  if(is.character(map) && length(map) == 1)
    map = loadMap(map)
  else if(isChromMap(map) || (is.list(map) && all(sapply(map, isChromMap))))
    map = genomeMap(map)
  else if(!isGenomeMap(map))
    stop2("Argument `map` must be either a `genomeMap`, a single `chromMap` or a single character")
  
  mapchrom = sapply(map, attr, "chrom")
  if(any(mapchrom == "X"))
    stop2("X chromosomal simulations are put on hold, but will be back in the near future.")
  
  # Model: Either "chi" or "haldane"
  model = match.arg(model)
  
  # Skip recombination
  if(is.null(skipRecomb) && !is.null(ids) && !setequal(ids, labels(x))) {
    FOU = founders(x)
    useids = if(length(ids) == 1) parents(x, ids) else ids
    fous = unlist(lapply(useids, function(id) intersect(FOU, ancestors(x, id, inclusive = TRUE))))
    fousUniq = unique.default(fous)
    
    # Skip founders who are ancestors of at most one ids
    counts = sapply(fousUniq, function(f) sum(fous == f))
    skipRecomb = fousUniq[counts == 1]
  }
    
  if (verbose) {
    message(glue::glue("
      Simulation parameters:
      Simulations  : {N}
      Chromosomes  : {toString2(mapchrom)}
      Genome length: {round(physRange(map), 2)} Mb
                     {round(mapLen(map, 'male'), 2)} cM (male)
                     {round(mapLen(map, 'female'), 2)} cM (female)
      Recomb model : {model}
      Target indivs: {toString2(ids, ifempty = '-')}
      Skip recomb  : {toString2(skipRecomb, ifempty = '-')}"
    ))
  }
  
  # Seed for random sampling
  if(!is.null(seed))
    set.seed(seed)
  
  # Start-data with founder alleles
  startData = distributeFounderAlleles(x, chrom = "AUTOSOMAL")
  
  # The actual simulations
  genomeSimList = replicate(N, 
    genedrop(x, ids, map = map, model = model, skipRecomb = skipRecomb, startData = startData),
    simplify = FALSE
  )
  
  # Timing
  if(verbose)
    message("Total time used: ", format(Sys.time() - starttime, digits = 3))
  
  structure(genomeSimList, 
            pedigree = x, 
            ids = ids, 
            skipRecomb = skipRecomb, 
            chrom = mapchrom, 
            physRange = physRange(map),
            mapLen = mapLen(map),
            model = model, 
            class = "genomeSimList")
}



#' @export
print.genomeSimList = function(x, ...) {
  attrs = attributes(x)
  len = attrs$mapLen
  
  print(glue::glue("
  List of {length(x)} genome simulations.
  Chromosomes: {toString2(attrs$chrom)}
  Total range: {round(attrs$physRange, 2)} Mb
  Map length : {len[[1]]} cM (male), {len[[2]]} cM (female)
  Rec. model : {attrs$model}
  Target ids : {toString2(attrs$ids)}
  Skip recomb: {toString2(attrs$skipRecomb, ifempty = '-')}
  "))
}

