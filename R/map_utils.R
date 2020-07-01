
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



#' Uniform recombination maps
#'
#' Create a uniform recombination map of a given length.
#'
#' @param Mb Map length in megabases.
#' @param cM Map length in centiMorgan.
#' @param M Map length in Morgan.
#' @param cm.per.mb A positive number; the cM/Mb ratio.
#' @param chrom A chromosome label.
#'
#' @return An object of class `chromMap`, which is a list of two matrices,
#'   named "male" and "female".
#'
#' @examples
#' uniformMap(M = 1)
#'
#' m = uniformMap(Mb = 1, cM = 2:3)
#' @export
uniformMap = function(Mb = NULL, cM = NULL, M = NULL, cm.per.mb = 1, 
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
    cM = if (!is.null(M)) M * 100 else cm.per.mb * Mb
  
  if (is.null(Mb)) 
    Mb = cM / cm.per.mb
  
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
#' This function loads one of the built-in genetic maps. A faster, uniform
#' version is also available by the parameter `detailed = FALSE`.
#'
#' @param map The name of the wanted map. Currently, the only valid choice is
#'   "decode19". This is also the default.
#' @param chrom A numeric vector indicating which chromosomes to load. Default:
#'   `1:22` (the autosomes)
#' @param detailed A logical. If TRUE (default), the complete inhomogeneous map
#'   is used. If FALSE, a uniform version of the same map is produced, i.e. with
#'   the correct lengths, but constant recombination rate along each chromosome.
#' @param sexSpecific A logical, by default TRUE. If FALSE, a sex-averaged map
#'   is returned, equal between males and females
#'
#' @return An object of class `genomeMap`.
#' 
#' @references Halldorsson et al. _Characterizing mutagenic effects of recombination through a
#'   sequence-level genetic map._ Science 363, no. 6425 (2019).
#' 
#' @examples
#' # By default, the complete map of all 22 autosomes is returned
#' loadMap()
#'
#' # Uniform version
#' m = loadMap(detailed = FALSE)
#'
#' # Check chromosome 1:
#' m1 = m[[1]]
#' m1$male
#' m1$female
#' 
#' @export
loadMap = function(map = "decode19", chrom = 1:22, detailed = TRUE, sexSpecific = TRUE) {
  
  if(!is.character(map) || length(map) != 1)
    stop2("Argument `map` must be a character of length 1")
  
  if(!is.logical(detailed) || length(detailed) != 1 || is.na(detailed))
    stop2("Argument `detailed` must be either TRUE or FALSE")

  if(!is.logical(sexSpecific) || length(sexSpecific) != 1 || is.na(sexSpecific))
    stop2("Argument `sexSpecific` must be either TRUE or FALSE")
  
  # For now only 'decode19' is implemented
  builtinMaps = c("decode19")
  mapno = pmatch(map, builtinMaps)
  if(is.na(mapno))
    stop2("Unknown map: ", map)

  map = builtinMaps[mapno]
  genome = get(map)[chrom]
  
  if(detailed && sexSpecific)
    return(genome)
  
  if(!detailed) {
    chroms = lapply(genome, function(chr) {
      chrom = attr(chr, "chrom")
      mb = chromLen(chr, "Mb")
      cm = chromLen(chr, "cM", sex = NA)
      if(!sexSpecific) 
        cm = mean(cm)
      
      uniformMap(Mb = mb, cM = cm, chrom = chrom)
    })
    
    return(genomeMap(chroms))
  }
  
  if(!sexSpecific) {
    genome[] = lapply(genome, function(chr) {
      if(!identical(chr$male$Mb, chr$female$Mb))
        stop2("Sex averaging requires equal map positions in males and females")
      
      chr$male$cM = chr$female$cM = (chr$male$cM + chr$female$cM)/2
      chr
    })
    
    return(genome)
  }
}


cm2phys = function(cM, map) { # map: columns 'Mb' and 'cM'
  if(!length(cM)) 
    return(cM)
  mapMB = map$Mb
  mapCM = map$cM
  
  res = numeric(length(cM))
  nontriv = cM >= 0 & cM <= mapCM[length(mapCM)]
  res[!nontriv] = NA
  
  cm = cM[nontriv]
  interv = findInterval(cm, mapCM, all.inside = TRUE)
  res[nontriv] = mapMB[interv] + (cm - mapCM[interv]) *
    (mapMB[interv + 1] - mapMB[interv]) / (mapCM[interv + 1] - mapCM[interv])
  res
}


phys2cm = function(Mb, map) {    # map: columns 'Mb' and 'cM'
  if(!length(Mb)) 
    return(Mb)
  
  mapMB = map$Mb
  mapCM = map$cM
  
  nontriv = Mb >= 0 & Mb <= mapMB[length(mapMB)]
  res = numeric(length(Mb))
  res[!nontriv] = NA
  
  mb = Mb[nontriv]
  interv = findInterval(mb, mapMB, all.inside = TRUE)
  res[nontriv] = mapCM[interv] + (mb - mapMB[interv]) *
    (mapCM[interv + 1] - mapCM[interv]) / (mapMB[interv + 1] - mapMB[interv])
  res
}
