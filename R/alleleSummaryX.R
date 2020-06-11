#' X-chromosomal allele sharing summary
#'
#' This is the X-chromosomal version of [alleleSummary()].
#'
#' @param x An object of class `genomeSim`, i.e. a list of simulated
#'   chromosomes. Each chromosome is a list, with one entry for each individual.
#'   Each of these entries is a list of two matrices (one for each strand). The
#'   matrices have 2 columns (start position; allele) and one row for each
#'   segment unbroken by recombination.
#' @param ids A vector of ID labels. If missing, all individuals are included.
#'
#' @return A numerical matrix. Each row corresponds to a chromosomal segment.
#'   The first 4 columns describe the segment (chromosome, start, end, length).,
#'   and are followed by one or two columns for each of the selected
#'   individuals: One column (maternal allele) for males and two columns
#'   (paternal allele, maternal allele) for females. If `length(ids) == 2` two
#'   additional columns are added:
#'
#'   * `IBD` : The IBD status of each segment (= number of alleles shared
#'   identical by descent). For a given segment, the IBD status is either 0, 1,
#'   2 or NA. If either individual is X-inbred, they may be autozygous in a
#'   segment, in which case the IBD status is reported as NA. With inbred
#'   individuals the `Sigma` column (see below) is more informative than the
#'   `IBD` column.
#'
#'   * `Sigma` : The condensed identity state of each segment, given as an
#'   integer in the range 1-9. The numbers refer to the _autosomal_ states in
#'   the usual ordering; for details about this, and how how they relate to
#'   identity states on X, please see the explanation at the `ribd` homepage:
#'   <https://github.com/magnusdv/ribd>.
#'
#' @examples
#' \donttest{
#' x = fullSibMating(1)
#' s = ibdsim(x, N = 1, ids = NULL, chrom = "X")[[1]]
#'
#' # Complete summary.
#' # Note only one allele column for males
#' alleleSummaryX(s)
#'
#' # Outbred brother/sister
#' alleleSummaryX(s, ids = 3:4)
#'
#' # Inbred brother/sister
#' alleleSummaryX(s, ids = 5:6)
#' }
#' @export
alleleSummaryX = function(x, ids) {
  
  ped = attr(x, "pedigree")
  if (missing(ids)) 
    ids = labels(ped)
  sex = getSex(ped, ids) 
  
  # First make the autosomal summary
  res = alleleSummary(x, ids)
  
  # Remove paternal columns of males
  male_pat = colnames(res) %in% paste0(ids[sex == 1], ":p")
  res = res[, !male_pat, drop = FALSE]
  
  # Recalculate IBD column
  if(length(ids) == 2 && sum(sex) < 4) {
    sigma = res[, "Sigma"]
    if(sex[1] == 1 && sex[2] == 1)
      IBD = 2 - sigma   # always either 1 or 2
    else if(sex[1] == 1 && sex[2] == 2)
      IBD = ifelse(sigma == 4, 0L, ifelse(sigma == 3, 1L, NA_integer_))
    else if(sex[1] == 2 && sex[2] == 1)
      IBD = ifelse(sigma == 6, 0L, ifelse(sigma == 5, 1L, NA_integer_))
    
    res[, 'IBD'] = IBD
  }
  res
}

