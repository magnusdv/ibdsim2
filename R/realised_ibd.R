#' Realised relatedness
#'
#' Compute the realised values of various measures of pairwise relatedness, in
#' simulated data. For now, only the realised IBD ("kappa") coefficients are implemented.
#'
#' DRAFT:  
#' Consider two members A and B of a pedigree P. The *kinship coefficient*
#' between A and B is defined as the probability that a random allele sampled in
#' A is identical by descent (IBD) with an allele sampled in B at the same
#' autosomal locus. If this probability is taken conditional only on the
#' pedigree P, the result is the traditional pedigre-based kinship coefficient.
#'
#' However, because of the discrete nature of meiotic recombination, the actual
#' IBD distribution of individuals with the specified relationship is subject to
#' variation. Hence we may also be interested in the *realised* (or
#' *genomic*) kinship coefficient between A and B, which is the same
#' probability as above, but conditional on the recombination events in the
#' meioses between A and B in P. For example, if the recombition events between
#' A and B happen to result in no segments of IBD sharing, the realised kinship
#' is 0, whatever the pedigree-based kinship may be.
#'
#' 
#' @param sim A list of genome simulations, as output by [ibdsim()].
#' @param id.pair A vector of length 2, with ID labels of the two individuals in question.
#'
#' @examples
#' x = pedtools::nuclearPed(2)
#' s = ibdsim(x, sims=10)
#' realised_kappa(s, id.pair=3:4)
#' 
#' @export
realised_kappa = function(sim, id.pair) {
  if(length(id.pair) != 2) 
    stop2("`id.pair` must be a vector of length 2")
  
  L = attr(sim, 'genome_length_Mb')
  
  segment_summary = vapply(sim, function(s) {
    a = alleleSummary(s, ids=id.pair, ibd.status=T)
    chrom = a[,'chrom']
    ibd = a[,'ibd']  

    # merge adjacent segments with equal IBD status (and equal chrom)
    seg_starts_idx = which(c(T, diff(ibd) != 0 | diff(chrom) != 0))
    seg_ends_idx = c(seg_starts_idx[-1] - 1, length(ibd))
    
    a_merged = a[seg_starts_idx, c('chrom', 'start', 'end', 'length', 'ibd')]
    a_merged[, 'end'] = a[seg_ends_idx, 'end']
    a_merged[, 'length'] = a_merged[, 'end'] - a_merged[, 'start'] 
    # TODO: Possible speedup of the above: Modify 'end' and 'length' only when needed
    
    len = a_merged[, 'length']  
    ibd = a_merged[, 'ibd']  
      
    c(ibd0 = sum(len[ibd==0]), ibd1 = sum(len[ibd==1]), ibd2 = sum(len[ibd==2]), 
      Nseg1 = sum(ibd==1), Nseg2 = sum(ibd==2), Nseg= sum(ibd>0))
  }, numeric(6))
  
  kappa.realised = segment_summary[1:3, , drop = FALSE]/L
  list(kappa.realised = kappa.realised, 
       Nsegments = segment_summary[4:6, , drop = FALSE], 
       kappa.hat= rowMeans(kappa.realised),
       genomeLength = L)
}
