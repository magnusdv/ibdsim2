#' Uniform recombination maps
#'
#' Create a uniform recombination map of a given length.
#'
#' @param Mb Map length in megabases.
#' @param cM Map length in centiMorgan.
#' @param M Map length in Morgan.
#' @param cmPerMb A positive number; the cM/Mb ratio.
#' @param chrom A chromosome label.
#'
#' @return An object of class `chromMap`, which is a list of two matrices,
#'   named "male" and "female".
#'
#' @seealso [loadMap()], [customMap()]
#' 
#' @examples
#' uniformMap(M = 1)
#'
#' m = uniformMap(Mb = 1, cM = 2:3)
#' @export
uniformMap = function(Mb = NULL, cM = NULL, M = NULL, cmPerMb = 1, 
                      chrom = 1) {
  
  if(is.null(cM) &&  is.null(M) && is.null(Mb))
    stop2("No map length indicated")
  
  if(!is.null(cM) && !is.null(M)) 
    stop2("Either `cM` or `M` must be NULL")
  
  if(!is.null(Mb) && !(is.numeric(Mb) && length(Mb) == 1))
    stop2("When non-NULL, `Mb` must be a numeric of length 1: ", Mb)
  
  if(!is.null(cM) && !(is.numeric(cM) && length(cM) < 3))
    stop2("When non-NULL, `cM` must be a numeric of length 1 or 2: ", cM)
  
  if(!is.null(M) && !(is.numeric(M) && length(M) < 3))
    stop2("When non-NULL, `M` must be a numeric of length 1 or 2: ", M)
  
  if (is.null(cM))
    cM = if (!is.null(M)) M * 100 else cmPerMb * Mb
  
  if (is.null(Mb)) 
    Mb = cM / cmPerMb
  
  # If length 0, return early
  if(Mb == 0) {
    male = female = cbind(Mb = 0, cM = 0)
    return(chromMap(male, female, chrom = chrom))
  }
  
  Mb = unname(rep(Mb, length.out = 2))
  cM = unname(rep(cM, length.out = 2))
  
  male = cbind(Mb = c(0, Mb[1]), cM = c(0, cM[1]))
  female = cbind(Mb = c(0, Mb[2]), cM = c(0, cM[2]))
  
  if (is.character(chrom) && tolower(chrom) == "x") {
    chrom = "X"
    male = NULL
  }
  
  chromMap(male, female, chrom = chrom)
}


#' Load a built-in genetic map
#'
#' This function loads one of the built-in genetic maps. Currently, the
#' available map is based on the publication by Halldorsson et al. (2019).
#'
#' For reasons of speed and efficiency, the built-in map is a thinned version of
#' the published map (Halldorsson et al., 2019), keeping around 60 000 data
#' points.
#'
#' By setting `uniform = TRUE`, a uniform version of the map is returned, in
#' which each chromosome has the same genetic lengths as in the original, but
#' with constant recombination rates. This gives much faster simulations and may
#' be preferable in some applications.
#'
#' @param map The name of the wanted map, possibly abbreviated. Currently, the
#'   only valid choice is "decode19" (default).
#' @param chrom A numeric vector indicating which chromosomes to load. Default:
#'   `1:22` (the autosomes).
#' @param uniform A logical. If FALSE (default), the complete inhomogeneous map
#'   is used. If TRUE, a uniform version of the same map is produced, i.e. with
#'   the correct lengths, but constant recombination rate along each chromosome.
#' @param sexAverage A logical, by default FALSE. If TRUE, a sex-averaged map is
#'   returned, with equal recombination rates for males and females.
#'
#' @return An object of class `genomeMap`.
#'
#' @seealso [uniformMap()], [customMap()]
#' 
#' @references Halldorsson et al. _Characterizing mutagenic effects of
#'   recombination through a sequence-level genetic map._ Science 363, no. 6425
#'   (2019).
#'
#' @examples
#' # By default, the complete map of all 22 autosomes is returned
#' loadMap()
#'
#' # Uniform version
#' m = loadMap(uniform = TRUE)
#'
#' # Check chromosome 1
#' m1 = m[[1]]
#' m1$male
#' m1$female
#'
#' @export
loadMap = function(map = "decode19", chrom = 1:22, uniform = FALSE, sexAverage = FALSE) {
  
  if(!is.character(map) || length(map) != 1)
    stop2("Argument `map` must be a character of length 1")
  
  if(!is.logical(uniform) || length(uniform) != 1 || is.na(uniform))
    stop2("Argument `uniform` must be either TRUE or FALSE")
  
  if(!is.logical(sexAverage) || length(sexAverage) != 1 || is.na(sexAverage))
    stop2("Argument `sexAverage` must be either TRUE or FALSE")
  
  # For now only 'decode19' is implemented
  builtinMaps = c("decode19")
  mapno = pmatch(map, builtinMaps)
  if(is.na(mapno))
    stop2("Unknown map: ", map)
  
  map = builtinMaps[mapno]
  genome = get(map)[chrom]
  
  if(!uniform && !sexAverage)
    return(genome)
  
  if(uniform) {
    chroms = lapply(genome, function(chr) {
      chrom = attr(chr, "chrom")
      mb = chromLen(chr, "Mb")
      cm = chromLen(chr, "cM", sex = NA)
      
      if(sexAverage) 
        cm = mean(cm)
      
      uniformMap(Mb = mb, cM = cm, chrom = chrom)
    })
    
    return(genomeMap(chroms))
  }
  
  if(sexAverage) {
    genome[] = lapply(genome, function(chr) {
      if(!identical(chr$male$Mb, chr$female$Mb))
        stop2("Sex averaging requires equal map positions in males and females")
      
      chr$male$cM = chr$female$cM = (chr$male$cM + chr$female$cM)/2
      chr
    })
    
    return(genome)
  }
}


#' Custom recombination map
#'
#' Create custom recombination maps for use in [ibdsim()].
#'
#' The column names of `x` must include either
#'
#' * `chrom`, `mb` and `cm`   (sex-averaged map)
#'
#' or
#'
#' * `chrom`, `mb`, `male` and `female`  (sex-specific map)
#'
#' Upper-case letters are allowed in these names. The `mb` column should contain
#' physical positions in megabases, while `cm`, `male`, `female` give the
#' corresponding genetic position in centiMorgans.
#'
#' @param x A data frame or matrix. See details for format specifications.
#'
#' @return An object of class `genomeMap`.
#'
#' @seealso [uniformMap()], [loadMap()]
#'
#' @examples
#' # A map including two chromosomes.
#' df1 = data.frame(chrom = c(1, 1, 2, 2),
#'                  mb = c(0, 2, 0, 5),
#'                  cm = c(0, 3, 0, 6))
#' map1 = customMap(df1)
#' map1
#'
#' # Use columns "male" and "female" to make sex specific maps
#' df2 = data.frame(chrom = c(1, 1, 2, 2),
#'                  mb = c(0, 2, 0, 5),
#'                  male = c(0, 3, 0, 6),
#'                  female = c(0, 4, 0, 7))
#' map2 = customMap(df2)
#' map2
#'
#' @export
customMap = function(x) {
  if(is.matrix(x))
    x = as.data.frame(x)
  
  if(!is.data.frame(x))
    stop2("Argument `x` must be a data frame or matrix. Received: ", class(x))
  
  nms = names(x)
  nmsSmall = tolower(nms)
  if(!"chrom" %in% nmsSmall) stop2('`x` must have a column named "chrom". See `?customMap`')
  if(!"mb" %in% nmsSmall) stop2('`x` must have a column named "mb". See `?customMap`')
  sexEq = "cm" %in% nmsSmall
  sexSpec = "male" %in% nmsSmall && "female" %in% nmsSmall
  if(!(sexEq || sexSpec))
    stop2('`x` must either have a colum named "cm", or two columns named "male" and "female". See `?customMap`')
  
  names(x) = nmsSmall
  
  # Split by chromosome
  chrList = split(x, x$chrom)
  
  chroms = lapply(chrList, function(chr) {
    if(sexEq) {
      male = female = chr[c("mb", "cm")]
    }
    else {
      male = chr[c("mb", "male")]
      female = chr[c("mb", "female")]
    }
    chromMap(male, female, chrom = chr$chrom[1])
  })
  
  genomeMap(chroms)
}


chromMap = function(male, female = male, chrom = 1) {
  dmm = dim(male)
  dmf = dim(female)
  if(is.null(dmm) || dmm[2] != 2)
    stop2("Male map does not have two columns")
  if(is.null(dmf) || dmf[2] != 2)
    stop2("Female map does not have two columns")
  
  # Convert matrix/tibbles/etc
  male = as.data.frame(male)
  female = as.data.frame(female)
  
  # Fix names
  names(male) = names(female) = c("Mb", "cM")
  
  physM = male$Mb
  physF = female$Mb
  
  if(physM[1] != physF[1])
    stop2("First position must be the same in male and female maps: ", c(physM[1], physF[1]))
  
  if(physM[dmm[1]] != physF[dmf[1]])
    stop2("End position must be the same in male and female maps: ", c(physM[dmm[1]], physF[dmf[1]]))
  
  if(physF[dmf[1]] > 1e9) {
    male$Mb = male$Mb / 1e6
    female$Mb = female$Mb / 1e6
  }
  
  physLen = female$Mb[dmf[1]] - female$Mb[1]
  
  structure(list(male = male, female = female), chromLen = physLen, chrom = chrom, 
            class = "chromMap")
}

genomeMap = function(x) {
  if(isChromMap(x))
    x = list(x)
  else if(!all(sapply(x, isChromMap)))
    stop2("Input to `genomeMap()` must a list of `chromMap' objects")
  
  len = sum(sapply(x, chromLen))
  structure(x, genomeLen = len, class = "genomeMap")
}


#' @export
`[.genomeMap` = function(x, i) {
  if(!all(i %in% seq_along(x)))
    stop2("Index out of range: ", setdiff(i, seq_along(x)))
  
  s = unclass(x)[i]
  len = sum(sapply(s, chromLen))
  structure(s, genomeLen = len, class = "genomeMap")
}

#' @export
print.chromMap = function(x, ...) {
  chr = attr(x, 'chrom')
  maLen = length(x$male$Mb)
  feLen = length(x$female$Mb)
  nPoints = if(maLen == feLen) feLen else sprintf("%d (male); %d (female)", maLen, feLen)
  
  print(glue::glue("
  Map of chromosome {chr}
  Physical length: {round(chromLen(x), 2)} Mb
  Physical range : {paste(round(range(x$female$Mb), 2), collapse = ' - ')} Mb 
  Male length    : {round(chromLen(x, 'cM', sex = 'male'), 2)} cM
  Female length  : {round(chromLen(x, 'cM', sex = 'female'), 2)} cM
  Data points    : {nPoints}
  "))
}

#' @export
print.genomeMap = function(x, ...) {
  nChr = length(x)
  plural = if(nChr > 1) 's' else ''
  
  print(glue::glue("
  Genome map consisting of {nChr} chromosome{plural}
  Physical length: {round(genomeLen(x), 2)} Mb
  Male length    : {round(genomeLen(x, 'cM', sex = 'male'), 2)} cM
  Female length  : {round(genomeLen(x, 'cM', sex = 'female'), 2)} cM
  "))
}
  
isChromMap = function(x)
  inherits(x, "chromMap")

isGenomeMap = function(x)
  inherits(x, "genomeMap")

chromLen = function(x, unit = c("Mb", "cM"), sex = NA) {
  if(!isChromMap(x))
    stop2("First argument must be a `chromMap` object, not ", class(x))
  
  unit = match.arg(unit)
  if(unit == "Mb")
    res = attr(x, "chromLen") %||% attr(x, "length_Mb")
  else {
    if(is.na(sex)) {
      ma = chromLen(x, 'cM', 'male')
      fe = chromLen(x, 'cM', 'female')
      return(c(male = ma, female = fe))
    }
    
    df = x[[match.arg(sex, c('male', 'female'))]]
    res = df$cM[nrow(df)]
  }
  if(is.null(res)) 
    res = 0
  
  res
}

genomeLen = function(x, unit = c("Mb", "cM"), sex = NA) { 
  if(!isGenomeMap(x))
    stop2("First argument must be a `genomeMap` object, not ", class(x))
  
  unit = match.arg(unit)
  if(unit == "Mb")
    res = attr(x, "genomeLen") %||% attr(x, "length_Mb") %||% attr(x, "genome_length_Mb")
  else {
    lens = sapply(x, chromLen, unit = "cM", sex = sex)
    if(is.na(sex))
      res = rowSums(lens)
    else 
      res = sum(lens)
  }
  
  res
}




#' Conversion of genetic map positions
#'
#' Convert between physical position (in megabases) and genetic position
#' (centiMorgan) given a chromosome map. Linear extrapolation is used to convert
#' positions between map points.
#'
#' @param Mb A vector of physical positions (in Mb), or NULL.
#' @param cM A vector of genetic positions (in cM), or NULL.
#' @param map A data frame with columns `Mb` and `cM`.
#'
#' @return A vector of the same length as the input.
#'
#' @examples
#' # Chromosome 1 of the built-in recombination map
#' map = loadMap(chrom = 1)[[1]]
#' head(map$male)
#'
#' # Conversion Mb -> cM
#' phys = 1:5
#' gen = convertPos(Mb = phys, map = map$male)
#' gen
#'
#' # Convert back (note the first position, which was outside of map)
#' convertPos(cM = gen, map = map$male)
#'
#' @export
convertPos = function(Mb = NULL, cM = NULL, map) {
  if(is.null(Mb) + is.null(cM) != 1)
    stop2("Exactly one of `Mb` and `cM` must be NULL")
  
  if(is.null(Mb)) {
    from = cM
    fromUnit = "cM"
    toUnit = "Mb"
  } else {
    from = Mb
    fromUnit = "Mb"
    toUnit = "cM"
  }
  
  # Output vector
  res = numeric(length(from))
  
  mapFrom = map[[fromUnit]]
  mapTo = map[[toUnit]]
  mapN = nrow(map)
  
  # Deal with values outside of map
  before = from < mapFrom[1]
  after =  from > mapFrom[mapN]
  res[before] = 0
  res[after] = if(toUnit == "Mb") NA else mapTo[mapN]
  
  # Locate inside values
  inside = !before & !after
  src = from[inside]
  idx = findInterval(src, mapFrom, all.inside = TRUE)
  
  # Rate of change
  rate = (mapTo[idx + 1] - mapTo[idx]) / (mapFrom[idx + 1] - mapFrom[idx])
  
  # Extrapolate
  res[inside] = mapTo[idx] + (src - mapFrom[idx]) * rate
  
  res
}
