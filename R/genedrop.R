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
  Xchrom = identical(chrom, "X")
  Xmale = Xchrom & x$SEX == 1
  
  maplen = attr(chrommap, "mapLen") # length in cM
  startMb = attr(chrommap, "physStart")
  endMb = attr(chrommap, "physEnd")
  
  skip = x$ID %in% skipRecomb
  
  h = startData %||% distributeFounderAlleles(x, Xchrom = Xchrom)
  
  for(i in NONFOU) {
    fa = FIDX[i]
    mo = MIDX[i]
    matGamete = meiosis(h[[mo]], map = chrommap$female, maplen = maplen[2], 
                        model = model, skipRecomb = skip[mo])
    if(Xmale[i])
      patGamete = NULL
    else
      patGamete = meiosis(h[[fa]], map = chrommap$male, maplen = maplen[1], 
                          model = model, skipRecomb = skip[fa])
    
    h[[i]] = list(pat = patGamete, mat = matGamete)
  }
  
  ### Convert to matrix of allele segments ###
  haplos = unlist(h[IDS], recursive = FALSE, use.names = FALSE)
  breaks = unlist(lapply(haplos, function(m) m[-1, 1]), use.names = FALSE)
  if(anyDuplicated.default(breaks))
    breaks = unique.default(breaks)
  sta = c(startMb, .sortDouble(breaks))
  
  # New, fast version
  alleleMat = build_allelemat_C(sta, haplos)
  
  # Previous, slower:
  # alleleMat = vapply(haplos, function(m) pos2allele(m, posvec = sta), FUN.VALUE = sta)
  # if (length(sta) == 1)      # since vapply simplifies if FUN.VALUE has length 1
  #   dim(alleleMat) = c(1, 2 * length(IDS))
  
  sto = c(sta[-1], endMb)
  cbind(chrom = if(Xchrom) 23L else chrom, start = sta, end = sto, length = sto - sta, alleleMat)
}



# Start-data for genedrop
# List of length pedsize, with entries to become list(pat, mat).
distributeFounderAlleles = function(x, Xchrom = FALSE) {
  
  fou = founders(x, internal = TRUE)
  nfou = length(fou)
  
  # Auxiliary function producing a full haplotype (1x2 matrix; pos-allele) of allele `a`
  TMP = rbind(c(0,0))
  .setAllele = function(a) `[<-`(TMP,2,a)
  
  # Logical of length nfou
  Xmale = Xchrom & x$SEX[fou] == 1
  
  # Logical of length nfou (100% inbred founders?)
  FOU_INB = founderInbreeding(x, chromType = if(Xchrom) "x" else "autosomal")
  inb1 = FOU_INB > 0
  if(any(FOU_INB[inb1] < 1))
      stop2("Founder inbreeding less than 100% is not supported: ", setdiff(FOU_INB, c(0,1)))
  
  # Create output list
  h = vector("list", pedsize(x))
  names(h) = x$ID
  
  # Fill in founders
  for(i in seq_along(fou)) {
    if(Xmale[i]) {
      pat = NULL
      mat = .setAllele(2*i)
    }
    else {
      pat = .setAllele(2*i - 1)
      mat = if(inb1[i]) pat else .setAllele(2*i)
    }
    h[[fou[i]]] = list(pat = pat, mat = mat)
  }
  h
}
