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
#' @return An object of class `chromosomeMap`, which is a list of two matrices,
#'   named "male" and "female".
#'
#' @examples
#' uniformMap(M = 1)
#'
#' @export
uniformMap = function(Mb = NULL, cM = NULL, M = NULL, cm.per.mb = 1, 
                      chrom = 1) {
  
  if (!is.null(cM) && !is.null(M)) 
    stop2("Either `cM` or `M` must be NULL")
  stopifnot(!is.null(cM) || !is.null(M) || !is.null(Mb))

  if (is.null(cM))
    cM = if (!is.null(M)) M * 100 else  cm.per.mb * Mb
  if (is.null(Mb)) Mb = cM / cm.per.mb

  if (is.character(chrom) && tolower(chrom) == "x")
    chrom = "X"

  map = switch(max(length(Mb), length(cM)), {
    m = cbind(Mb = c(0, Mb), cM = c(0, cM))
    list(male = m, female = m)
  }, {
    if (length(cM) == 1) cM = c(cM, cM)
    if (length(Mb) == 1) Mb = c(Mb, Mb)
    list(male = cbind(Mb = c(0, Mb[1]), cM = c(0, cM[1])), 
         female = cbind(Mb = c(0, Mb[2]), cM = c(0, cM[2])))
  })
  
  female_phys = as.numeric(map$female[2, 1])
  if (identical(chrom, "X"))
    map$male = NA
  else if (female_phys != map$male[2, 1]) 
    stop2("Male and female chromosomes must have equal physical length")
  
  structure(map, length_Mb = female_phys, chrom = chrom, class = "chromosomeMap")
}




chromMap = function(male, female = male, chrom = 1) {
  dmm = dim(male)
  dmf = dim(female)
  if(is.null(dmm) || dmm[2] != 2)
    stop2("Male map does not have two columns")
  if(is.null(dmf) || dmf[2] != 2)
    stop2("Female map does not have two columns")
  if(!is.numeric(male))
    stop2("Male map is not numeric")
  if(!is.numeric(female))
    stop2("Female map is not numeric")

  if(is.matrix(male)) 
    male = as.data.frame(male)
  if(is.matrix(female)) 
    female = as.data.frame(female)
  
  names(male) = names(female) = c("Mb", "cM")
  
  physM = male$Mb
  physF = female$Mb
  
  if(physM[1] != physF[1])
    stop2("First position must be the same in male and female maps: ", physM[1], physF[1])
  
  if(physM[dmm[1]] != physF[dmf[1]])
    stop2("End position must be the same in male and female maps: ", physM[dmm[1]], physF[dmf[1]])
  
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

isChromMap = function(x)
  inherits(x, "chromMap")

isGenomeMap = function(x)
  inherits(x, "genomeMap")

chromLen = function(x) 
  attr(x, "chromLen") %||% attr(x, "length_Mb")

genomeLen = function(x) 
  attr(x, "genomeLen") %||% attr(x, "genome_length_Mb")


loadMap = function(map, chrom = NULL) {

  if (is.character(map)) {
    if(!all(chrom %in% c(1:23, "X")))
      stop2("Chromosome not found in the Decode map: ", setdiff(chrom, c(1:23, "X")))
    
    CHROM.LENGTH = cbind(male_morgan = c(1.9, 1.752, 1.512, 1.35, 1.302, 1.162, 1.238, 1.089, 1.047, 1.147, 0.992, 1.154, 0.919, 0.857, 0.825, 0.88, 0.863, 0.737, 0.708, 0.563, 0.426, 0.45, NA),
    female_morgan = c(3.34, 3.129, 2.687, 2.6, 2.49, 2.333, 2.21, 2.042, 1.879, 2.062, 1.886, 1.974, 1.468, 1.24, 1.393, 1.494, 1.529, 1.379, 1.152, 1.109, 0.67, 0.662, 1.733),
    Mb = c(247.2, 242.7, 199.3, 191.1, 180.6, 170.8, 158.7, 146.2, 140.1, 135.3, 134.4, 132.3, 114.1, 105.3, 100.2, 88.7, 78.6, 76.1, 63.8, 62.4, 46.9, 49.5, 154.6))

    if (is.null(chrom) || identical(chrom, "AUTOSOMAL")) 
      chromnum = 1:22
    else if (identical(chrom, "X")) 
      chromnum = 23
    else if(is.numeric(chrom))
      chromnum = chrom

    maps = switch(tolower(map),
      decode = DecodeMap[chromnum],
      uniform.sex.spec = lapply(chromnum, function(chr) {
        dat = as.numeric(CHROM.LENGTH[chr, ])
        uniformMap(M = dat[1:2], Mb = dat[3], chrom = chr)
      }),
      uniform.sex.aver = lapply(chromnum, function(chr) {
        dat = as.numeric(CHROM.LENGTH[chr, ])
        uniformMap(M = mean(dat[1:2]), Mb = dat[3], chrom = chr)
      }),
      stop2("Invalid map name"))
    
    # Fix chrom attributes TODO: update DecodeMap
    maps = lapply(maps, function(m) {
      attrs = attributes(m)
      chr = attrs$chromosome %||% attrs$chrom
      attrs$chromosome = NULL
      attrs$chrom = if(chr == 23) "X" else chr
      attributes(m) = attrs
      m
    })
    
  }
  else {
    maps = map
    if (inherits(maps, "chromosomeMap")) maps = list(maps)
  }
  attr(maps, "length_Mb") = sum(unlist(lapply(maps, attr, "length_Mb")))
  maps
}


cm2phys = function(cM_locus, mapmat) { # mapmat matrise med kolonner 'Mb' og 'cM'
  if(!length(cM_locus)) 
    return(cM_locus)
  mapMB = mapmat[, 'Mb']
  mapCM = mapmat[, 'cM']
  
  res = numeric(length(cM_locus))
  nontriv = cM_locus >= 0 & cM_locus <= mapCM[length(mapCM)]
  res[!nontriv] <- NA
  
  cm = cM_locus[nontriv]
  interv = findInterval(cm, mapCM, all.inside = TRUE)
  res[nontriv] = mapMB[interv] + (cm - mapCM[interv]) *
    (mapMB[interv + 1] - mapMB[interv]) / (mapCM[interv + 1] - mapCM[interv])
  res
}


phys2cm = function(Mb, mapmat) {    # mapmat matrise med kolonner 'Mb' og 'cM'
  if(!length(Mb)) 
    return(Mb)
  
  mapMB = mapmat$Mb
  mapCM = mapmat$cM
  
  nontriv = Mb >= 0 & Mb <= mapMB[length(mapMB)]
  res = numeric(length(Mb))
  res[!nontriv] <- NA
  mb = Mb[nontriv]
  interv = findInterval(mb, mapMB, all.inside = TRUE)
  res[nontriv] = mapCM[interv] + 
    (mapCM[interv + 1] - mapCM[interv]) * (mb - mapMB[interv]) / (mapMB[interv + 1] - mapMB[interv])
  res
}
