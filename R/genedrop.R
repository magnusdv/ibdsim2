#' @importFrom stats runif
genedrop = function(x, map, model = "chi", skipRecomb = NULL) {
  FIDX = x$FIDX
  MIDX = x$MIDX
  FOU = founders(x, internal = TRUE)
  NONFOU = nonfounders(x, internal = TRUE)
  chrom = attr(map, "chrom")
  
  h = distributeFounderAlleles(x, chrom)
  
  if (chrom == 23) {
    for (i in NONFOU) {
      fa = FIDX[i]
      mo = MIDX[i]
      maternal.gamete = meiosis(h[[mo]], map = map$female, model = model, skipRecomb = mo %in% skipRecomb)
      paternal.gamete = if (x$SEX[i] == 1) maternal.gamete else h[[fa]][[1]]
      h[[i]] = list(paternal.gamene, maternal.gamete)
    }
  }
  else {
    for (i in NONFOU) {
      fa = FIDX[i]
      mo = MIDX[i]
      h[[i]] = list(meiosis(h[[fa]], map = map$male, model = model, skipRecomb = fa %in% skipRecomb),
                    meiosis(h[[mo]], map = map$female, model = model, skipRecomb = mo %in% skipRecomb)
      )
    }
  }
  
  attr(h, "chrom") = chrom
  attr(h, "length_Mb") = attr(map, "length_Mb")
  attr(h, "model") = model
  attr(h, "skipped") = skipRecomb
  
  class(h) = "chromosomeSim"
  h
}


distributeFounderAlleles = function(x, chrom = "AUTOSOMAL") {
  h = vector("list", pedsize(x))
  fou = founders(x, internal = TRUE)
  nfou = length(fou)
  
  if(identical(chrom, "X")) {
    SEX = x$SEX
    alleles = numeric(nfou)
    k = 1
    for (i in seq_along(fou)) {
      sex = SEX[fou[i]]
      alleles[c(2 * i - 1, 2 * i)] = c(k, k + sex - 1)
      k = k + sex
    }
    aux = cbind(rep.int(0, 2 * nfou), alleles, deparse.level = 0)
  }
  else {
    aux = cbind(rep.int(0, 2 * nfou), seq_len(2 * nfou))
    
    # For 100% inbred founders; make alleles identical
    FOU_INB = founderInbreeding(x)
    stopifnot(all(FOU_INB %in% c(0,1)))
    inb1 = which(FOU_INB == 1)
    if(length(inb1) > 0)
      aux[2 * inb1, ] = aux[2 * inb1 - 1, ]
  }
  
  h[fou] = lapply(2 * seq_along(fou), function(i) 
    list(aux[i - 1, , drop = FALSE], aux[i, , drop = FALSE]))
  h
}

