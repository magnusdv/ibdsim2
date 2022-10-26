genedrop = function(x, ids, map, model, skipRecomb, startData) {
  
  # Loop over chromosomes
  sims = lapply(map, function(chrommap)
    genedrop.singlechrom(x, ids, chrommap, model, skipRecomb, startData))
  
  res = do.call(rbind, sims)
  
  # Name allele columns (safely)
  acols = paste(rep(ids, each = 2), c("p", "m"), sep = ":") 
  colnames(res)[seq.int(to = ncol(res), along.with = acols)] = acols
  
  if(length(ids) < 3)
    res = addStates(res, acols = acols)
  
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
  
  physStart = attr(chrommap, "physStart")
  physEnd = attr(chrommap, "physEnd")
  
  sta = c(physStart, .sortDouble(breaks))
  sto = c(sta[-1], physEnd)
  
  # Allele matrix (fast version)
  alleleMat = build_allelemat_C(sta, haplos)
  
  # Previous, slower:
  # alleleMat = vapply(haplos, function(m) pos2allele(m, posvec = sta), FUN.VALUE = sta)
  # if (length(sta) == 1)      # since vapply simplifies if FUN.VALUE has length 1
  #   dim(alleleMat) = c(1, 2 * length(IDS))
  
  # Add columns with average CM positions
  avmap = chrommap$female
  if(!Xchrom) 
    avmap$cM = (chrommap$male$cM + chrommap$female$cM)/2
  startCM = .convertPos1(Mb = sta, map = avmap)
  endCM = .convertPos1(Mb = sto, map = avmap)
  
  # Collect as matrix
  cbind(chrom = if(Xchrom) 23L else chrom, 
        startMB = sta, endMB = sto,
        startCM = startCM, endCM = endCM,
        alleleMat)
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
