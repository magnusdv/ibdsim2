
#' @export
print.chromosomeSim = function(x, ...) {
  chrom = attr(x, 'chromosome')
  len = attr(x, 'length_Mb')
  model = attr(x, 'model')
  pedsize = length(x)
  skipped = attr(x, 'skipped')
  skip_str = if(!length(skipped)) "-" else toString(skipped)
  
  print(glue::glue("
  Simulation of a single chromosome
  Chrom name    : {chrom}
  Chrom length  : {len} Mb
  Recomb model  : {model}
  Pedigree size : {pedsize}
  Skipped recomb: {skip_str}
  "))
}

#' @export
print.genomeSim = function(x, ...) {
  attrs = attributes(x)
  if(!length(attrs$skipped)) attrs$skipped = "-"
  
  print(glue::glue("
  Total map length: {attrs$genome_length_Mb} Mb
  Chromosomes: {paste(attrs$chromosomes, collapse = ',')}
  Recombination model: {attrs$model}
  Pedigree members: {pedsize(attrs$ped)}
  Skipped recomb: {paste(attrs$skipped, collapse = ',')}
  "))
}

#' @export
print.genomeSimList = function(x, ...) {
  print(glue::glue("List of {length(x)} genome simulations."))
}
