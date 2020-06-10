#' IBD simulation
#'
#' This is the main function of the package, simulating the recombination
#' process in each meioses of a pedigree. The output summarises the IBD segments
#' between all or a subset of individuals.
#'
#' Each simulation starts by unique alleles being distributed to the pedigree
#' founders. In each meiosis, homologue chromosomes are made to recombine
#' according to a renewal process along the four-strand bundle, with chi square
#' distributed waiting times. (For comparison purposes, Haldane's Poisson model
#' for recombination is also implemented.)
#'
#' Recombination rates are sex-dependent, and vary along each chromosome
#' according to the recombination map specified by the `map` parameter. By
#' default, the complete Decode map of the human autosome is used (see
#' References). If `map = "uniform.sex.spec"`, the genetic chromosome *lengths*
#' are as in the Decode map, but the recombination rate is kept constant along
#' each chromosome. If `map = "uniform.sex.aver"`, sex averaged genetic
#' chromosome lengths are used (and constant recombination rates along each
#' chromosome).
#'
#' @param x A [pedtools::ped()] object.
#' @param N A positive integer indicating the number of simulations.
#' @param ids A subset of pedigree members whose IBD sharing should be analysed.
#'   If NULL, the simulations are returned unprocessed.
#' @param map The genetic map to be used in the simulations: One of the
#'   character strings "decode", "uniform.sex.spec", "uniform.sex.aver". (See
#'   Details.)
#' @param chrom A numeric vector indicating chromosome numbers, or either
#'   of the words "AUTOSOMAL" or "X". The default is 1:22, i.e., the human
#'   autosomes.
#' @param model Either "chi" (default) or "haldane", indicating the statistical
#'   model for recombination. (See details.)
#' @param skipRecomb An optional vector of ID labels indicating individuals
#'   whose meioses should be simulated without recombination. (Each child will
#'   then receive a random strand of each chromosome.) If NULL (default), nobody
#'   is skipped. Note that `skipRecomb = founders(x)` is a sensible (and time
#'   saving) choice in many applications.
#' @param seed An integer to be passed on to [set.seed()]).
#' @param verbose A logical.
#'
#' @return A list of `genomeSim` objects.
#'
#'   A `genomeSim` object is essentially a numerical matrix describing the
#'   allele flow through the pedigree in a single simulated. Each row
#'   corresponds to a chromosomal segment. The first 4 columns describe the
#'   segment (chromosome, start, end, length), and are followed by two columns
#'   (paternal allele, maternal allele) for each of the selected individuals. If
#'   `length(ids) == 2` two additional columns are added:
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
ibdsim = function(x, N = 1, ids = labels(x), map = "decode", chrom = NULL,
                  model = "chi", skipRecomb = NULL, 
                  seed = NULL, verbose = TRUE) {
  # Check input
  if(!is.ped(x))
    stop2("The first argument must be a `ped` object")
  if(!is_count(N))
    stop2("`N` must be a positive integer")
  if(!model %in% c("chi", "haldane"))
    stop2('Argument `model`` must be either "chi" or "haldane"')
  if(!all(founderInbreeding(x) %in% c(0,1)))
    stop2("Founder inbreeding coefficients other than 0 and 1 are not allowed")
  
  model_string = if(model == "chi") "Chi^2 renewal process" else "Haldane's poisson process"    
  
  # Ensure that parents precede their children
  if (!hasParentsBeforeChildren(x)) {
    if(verbose) message("Reordering so that all parents precede their children")
    x = parentsBeforeChildren(x)
  }
  
  # Start timer
  starttime = proc.time()

  # Load map and extract chromosome names.
  map = loadMap(map, chrom = chrom)
  mapchrom = attr(map, "chromosome") %||% sapply(map, attr, "chromosome")

  if (verbose) {
    skip_str = if (is.null(skipRecomb)) "-" else toString(skipRecomb)
    ids_str = if (is.null(ids)) "-" else toString(ids)
    
    message(glue::glue("
    No. of sims: {N}
    Chromosomes: {toString(mapchrom)}
    Rec. model : {model_string}
    Target ids : {ids_str}
    Skip recomb: {skip_str}
    "))
  }
  
  # Seed for random sampling
  set.seed(seed)
  
  # Various attributes which will be attached to the sims
  attribs = list(pedigree = x, 
                 ids = ids,
                 skipped = skipRecomb,
                 genome_length_Mb = attr(map, "length_Mb"),
                 chrom = mapchrom,
                 model = model_string)
  
  # The actual simulations: One sim at the time; each chromosome in turn 
  genomeSimList = lapply(1:N, function(i) {
    s = lapply(map, function(m)
      genedrop(x, map = m, model = model, skipRecomb = skipRecomb))
    attributes(s) = c(attribs, class = "genomeSim") 
    alleleSummary(s, ids)
  })
  
  if (verbose) {
    elapsed = (proc.time() - starttime)[["elapsed"]]
    message("Simulation finished in ", elapsed, " seconds.")
  }
  
  # Add attributes and class to the entire list 
  attributes(genomeSimList) = c(attribs, class = "genomeSimList")
  
  genomeSimList
}

