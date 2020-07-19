

#' @export
`[.genomeMap` = function(x, i) {
  if(!all(i %in% seq_along(x)))
    stop2("Index out of range: ", setdiff(i, seq_along(x)))
  
  s = unclass(x)[i]
  structure(s, class = "genomeMap")
}

#' @export
print.chromMap = function(x, ...) {
  attrs = attributes(x)
  maN = length(x$male$Mb)
  feN = length(x$female$Mb)
  nPoints = if(maN == feN) feN else sprintf("%d (male); %d (female)", maN, feN)
  
  print(glue::glue("
  Map of chromosome {attrs$chrom}
  Physical range: {round(attrs$physRange, 2)} Mb ({round(attrs$physStart, 2)} - {round(attrs$physEnd, 2)})
  Male length   : {round(mapLen(x, sex = 'male'), 2)} cM
  Female length : {round(mapLen(x, sex = 'female'), 2)} cM
  Data points   : {nPoints}
  "))
}

#' @export
print.genomeMap = function(x, ...) {
  nChr = length(x)
  plural = if(nChr > 1) 's' else ''
  
  print(glue::glue("
  Genome map consisting of {nChr} chromosome{plural}
  Physical range: {round(physRange(x), 2)} Mb
  Male length   : {round(mapLen(x, sex = 'male'), 2)} cM
  Female length : {round(mapLen(x, sex = 'female'), 2)} cM
  "))
}
  
isChromMap = function(x)
  inherits(x, "chromMap")

isGenomeMap = function(x)
  inherits(x, "genomeMap")


#' Physical and genetic map lengths
#'
#' Utility functions for extracting the physical or genetic length of chromosome
#' maps and genome maps.
#'
#' @param x A `chromMap` or `genomeMap` object.
#' @param sex Either "male", "female" or both.
#' @param ... Not used.
#'
#' @return `mapLen()` returns a numeric of the same length as `sex`, with the
#'   genetic length(s) in centiMorgan.
#'
#'   `physRange()` returns the physical length (in Mb) of the chromosome/genome
#'   covered by the map. For example, for a chromosome map starting at 2 Mb and
#'   ending at 8 Mb, the output is 6.
#' 
#' @seealso [loadMap()], [uniformMap()]
#' 
#' @examples
#' m = loadMap(chrom = 1:2)
#' m
#'
#' # Applied to `genomeMap` object:
#' physRange(m)
#' mapLen(m)
#'
#' # Applied to `chromMap` object:
#' physRange(m[[1]])
#' mapLen(m[[1]])
#'
#' @name maplengths
NULL

#' @rdname maplengths 
#' @export
mapLen = function(x, ...) {
  UseMethod("mapLen")
}

#' @rdname maplengths 
#' @export
mapLen.chromMap = function(x, sex = c("male", "female"), ...) {
  len = attr(x, "mapLen")
  if(length(sex) == 1)
    len[[sex]]
  else
    len[sex]
}

#' @rdname maplengths 
#' @export
mapLen.genomeMap = function(x, sex = c("male", "female"), ...) { 
  if(length(sex) == 1)
    sum(vapply(x, mapLen, sex = sex, FUN.VALUE = numeric(1)))
  else
    rowSums(vapply(x, mapLen, sex = sex, FUN.VALUE = numeric(2)))
}

#' @rdname maplengths 
#' @export
physRange = function(x, ...) {
  UseMethod("physRange")
}

#' @rdname maplengths 
#' @export
physRange.chromMap = function(x, ...) {
  attr(x, "physRange")
}

#' @rdname maplengths 
#' @export
physRange.genomeMap = function(x, ...) { 
  sum(vapply(x, physRange.chromMap, FUN.VALUE = numeric(1)))
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
