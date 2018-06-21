
#' @export
print.chromosomeSim = function(x, ...) {
  chrom = attr(x, 'chromosome')
  len = attr(x, 'length_Mb')
  model = attr(x, 'model')
  pedsize = length(x)
  skipped = attr(x, 'skipped')
  skip_str = if(!length(skipped)) "None" else paste(skipped, collapse=",")
  cond_sap = attr(x, 'condition')
  cond = if(is.null(cond_sap)) "No" else cond_sap
  
  print(glue::glue("
  Simulation of a single chromosome throughout a pedigree.
  Chromosome name: {chrom}
  Chromosome length: {len} Mb
  Recombination model: {model}
  Pedigree members: {pedsize}
  Skipped recombination in: {skip_str}
  Conditional: {cond}
  "))
}

#' @export
print.genomeSim = function(x, ...) {
  attrs = attributes(x)
  if(!length(attrs$skipped)) attrs$skipped = "None"
  if(is.null(attrs$condition)) attrs$condition = "No"
  
  print(glue::glue("
  Total map length: {attrs$genome_length_Mb} Mb
  Chromosomes: {paste(attrs$chromosomes, collapse=',')}
  Recombination model: {attrs$model}
  Pedigree members: {attrs$ped$NIND}
  Skipped recombination in: {paste(attrs$skipped, collapse=',')}
  Conditional: {attrs$condition}
  "))
}

#' @export
print.genomeSimList = function(x, ...) {
  print(glue::glue("List of {length(x)} genome simulations."))
  print(x[[1]])
}
