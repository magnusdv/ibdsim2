genedrop = function(x, ids, map, model, skipRecomb, startData) {
  
  # Loop over chromosomes
  sims = lapply(map, function(chrommap)
    genedrop.singlechrom(x, ids, chrommap, model, skipRecomb, startData))
  
  res = do.call(rbind, sims)
  colnames(res)[-(1:4)] = paste(rep(ids, each = 2), c("p", "m"), sep = ":")
  
  if(length(ids) < 3)
    res = addStates(res)
  
  structure(res, class = "genomeSim")
}

genedrop.singlechrom = function(x, ids, chrommap, model, skipRecomb, startData) {
  FIDX = x$FIDX
  MIDX = x$MIDX
  IDS = internalID(x, ids)
  NONFOU = nonfounders(x, internal = TRUE)
  
  chrom = attr(chrommap, "chrom")
  skip = x$ID %in% skipRecomb
  
  h = startData %||% distributeFounderAlleles(x, chrom)
  
  for (i in NONFOU) {
    fa = FIDX[i]
    mo = MIDX[i]
    matGamete = meiosis(h[[mo]], map = chrommap$female, model = model, skipRecomb = skip[mo])
    if (chrom == "X")
      patGamete = if (x$SEX[i] == 1) matGamete else h[[fa]][[1]]
    else
      patGamete = meiosis(h[[fa]], map = chrommap$male, model = model, skipRecomb = skip[fa])
    
    h[[i]] = list(patGamete, matGamete)
  }
  
  ### Convert to matrix of allele segments ###
  haplos = unlist(h[IDS], recursive = FALSE)
  breaks = unlist(lapply(haplos, function(m) m[-1, 1]))
  if(anyDuplicated.default(breaks))
    breaks = unique.default(breaks)
  
  sta = c(0, .sortDouble(breaks))
  
  alleleMat = vapply(haplos, pos2allele, posvec = sta, FUN.VALUE = sta)
  if (length(sta) == 1)      # since vapply simplifies if FUN.VALUE has length 1
    dim(alleleMat) = c(1, 2 * length(IDS))
  
  sto = c(sta[-1], chromLen(chrommap))
  cbind(chrom = chrom, start = sta, end = sto, length = sto - sta, alleleMat)
}


# Start-data for genedrop: founder alleles
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
