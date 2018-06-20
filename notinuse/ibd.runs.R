sap.segments = function(twostrands_list, sap) { # genedrops for several chromosomes. sap = single allele pattern
  do.call(rbind, lapply(twostrands_list, sap.finder, sap = sap))
}

sap.finder = function(x, sap) { # x a list of haplo-objects (i.e. each element is a list of two matrices); sap = single allele pattern
  zero = sap[["0"]]; one = sap[["1"]]; two = sap[["2"]]; atm1 = sap[["atmost1"]]; atl1 = sap[["atleast1"]]; not1 = sap[["not1"]]
  stopifnot(!is.null(id1 <- c(two, atl1, one)[1]))
  sap.indivs = c(zero, one, two, atm1, atl1, not1)
  breaks = unlist(lapply(unlist(x[sap.indivs], recursive = FALSE), function(m) m[-1, 1]))
  breaks = c(0, .sortDouble(breaks[!duplicated(breaks)]))
  running_alleles = numeric(0); res = matrix(nrow = 0, ncol = 3); starts = numeric();
  ALLalleles = list(); ALLalleles[sap.indivs] <- lapply(x[sap.indivs], .getAlleles, posvec = breaks)
  for (i in seq_along(breaks)) {
    a = ALLalleles[[id1]][1, i]; b = ALLalleles[[id1]][2, i]
    for (A in running_alleles[!running_alleles %in% c(a, b)]) { # reached end of segment, write to result matrix
      res = rbind(res, c(starts[A], breaks[i], A))
      running_alleles = running_alleles[running_alleles != A]
    }
    for (al in {
      if (a == b) a else c(a, b)
    }) {
      ibd_nr = numeric(length(x))
      ibd_nr[sap.indivs] = unlist(lapply(ALLalleles[sap.indivs], function(als) sum(als[, i] == al)))
      test = all(ibd_nr[two] == 2) && all(ibd_nr[atl1] > 0) && all(ibd_nr[one] == 1) && all(ibd_nr[atm1] < 2) && all(ibd_nr[zero] == 0) && all(ibd_nr[not1] != 1)  # longer but faster

      ind = match(al, running_alleles, nomatch = 0) #
      if (ind > 0 && !test) { # cat(al,"pos",breaks[i],"stop(fail)\n")
        res = rbind(res, c(starts[al], breaks[i], al))
        running_alleles = running_alleles[-ind] # cat(running_alleles, "running\n")
      }
      else if (ind == 0 && test) { # cat(al,"pos",breaks[i],"start\n")
        starts[al] = breaks[i];
        running_alleles = c(running_alleles, al)
      }
    }
  }
  for (A in running_alleles) res = rbind(res, c(starts[A], attr(x, "length_Mb"), A)) # last entries

  if (!is.null(disloc <- attr(x, "dis.locus")))
    disreg = as.numeric(res[, 1] <= attr(x, "dis.locus") & res[, 2] > attr(x, "dis.locus") & res[, 3] == attr(x, "dis.allel"))
  else disreg = rep.int(0, nrow(res))

  chrom = rep.int(attr(x, "chromosome"), nrow(res))
  res = cbind(chrom, res, disreg)
  colnames(res) = c("chrom", "start", "end", "allele", "disreg")
  res
}

### not used!
# ibd.finder.ordered <-
# function(x, ibd.ordered, nonibd.ordered) {
# # x: single sim of a single chrom, all indivs. A list of haplo-objects (i.e. each element is a list of two matrices)
# # ibd: vector of signed ID integers. +/- = paternal/maternal.
# # Example: ibd = c(3,-4) ---> (paternal of ind 3 == maternal of ind 4 signed)
# all.conditions = c(ibd.ordered, nonibd.ordered)
# stopifnot(all(sapply(all.conditions, length) >= 2))
# ids = sort.int(unique.default(abs(unlist(all.conditions))))
# haplos = unlist(x[ibd], recursive=FALSE)
# breaks = unlist(lapply(haplos, function(m) m[-1,1]))
# breaks = c(0, .sortDouble(breaks[!duplicated(breaks)]))
# alleles = do.call(cbind, lapply(haplos, pos2allele, posvec=breaks))
# colnames(alleles) = paste0(c("","-"), rep(ids, each=2))

# test = rep(T, length(breaks))
# for(vec in ibd.ordered) {
# this_alleles = alleles[, as.character(vec)]
# test = test & rowSums(this_alleles - this_alleles[, 1]) == 0
# }
# for(vec in nonibd.ordered) {
# this_alleles = alleles[, as.character(vec)]
# test = test & rowSums(abs(this_alleles - this_alleles[, 1])) > 0
# }
# stops = c(breaks[-1], attr(x, 'length_Mb'))
# chrom = rep.int(attr(x, 'chromosome'), sum(test))
# disreg = rep.int(0, sum(test)) # Dont want/need this? To make sure downstream summary stuff works...
# cbind(chrom=chrom, start=breaks[test], end=stops[test], alleles[test, , drop=F], disreg=disreg)
# }

### not used!!!
# ibd.status <-
# function(x, id.pair) {
# # x: single sim of a single chrom, all indivs. A list of haplo-objects (i.e. each element is a list of two matrices)
# haplos = unlist(x[id.pair], recursive=FALSE)
# breaks = unlist(lapply(haplos, function(m) m[-1,1]))
# breaks = c(0, .sortDouble(breaks[!duplicated(breaks)]))
# alleles_list = lapply(haplos, pos2allele, posvec=breaks)

# ibd_pp = alleles_list[[1]] == alleles_list[[3]]
# ibd_pm = alleles_list[[1]] == alleles_list[[4]]
# ibd_mp = alleles_list[[2]] == alleles_list[[3]]
# ibd_mm = alleles_list[[2]] == alleles_list[[4]]
# ibd = (ibd_pp | ibd_pm) + (ibd_mp | ibd_mm)

# starts = breaks
# stops = c(breaks[-1], attr(x, 'length_Mb'))
# alleles_matrix = do.call(cbind, alleles_list)
# colnames(alleles_matrix) = paste0(rep(id.pair, each=2), c("p","m"))
# chrom = rep.int(attr(x, 'chromosome'), length(breaks))

# cbind(chrom=chrom, start=starts, end=stops, length=stops-starts, alleles_matrix, ibd=ibd, ibd_pp=ibd_pp, ibd_pm=ibd_pm, ibd_mp=ibd_mp, ibd_mm=ibd_mm)
# }
