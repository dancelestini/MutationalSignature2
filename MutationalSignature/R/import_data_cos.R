#' Load a Cosmic targeted screens mutation data
#'
#' This function will load and delete all the no needed patient with NA values or Null characters.
#'
#' @param file Path to the input file
#' @return A dataframe with all usable patient
#' @export
import_data_cos = function(file){
  cosmic = read.csv(file,header = T)
  cosmic = cosmic[cosmic$MUTATION_AA != "null",]
  cosmic = cosmic[-which(cosmic$MUTATION_AA == "p.?"),]

  genes_mut_site1 = sapply(strsplit(cosmic$MUTATION_CDS,"[\\c.]+"), `[`, 2)
  change_nu = sapply(strsplit(genes_mut_site1,"[0-9]+"),tail,1)
  idx_del = which(change_nu %in% c("A>T","A>C","A>G","C>A","C>G","C>T","G>A","G>C","G>T","T>A","T>C","T>G"))
  cosmic = cosmic[idx_del,]
  return(cosmic)
}
