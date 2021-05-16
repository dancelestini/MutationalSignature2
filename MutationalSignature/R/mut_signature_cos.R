#' Create the final results Dataframe
#'
#' This function will create a dataframe with all possible mutation and count and calculate the frequency
#' of each from the input dataset.
#'
#' @param mut_to Codon of mutation paste to the change of nucleotide due to mutation.
#' @param file Just the number of the file analyzed. Default = 1
#' @return mut_signa = dataset with the counts of each mutations
#' @export
mut_signature_cos = function(mut_to,file = 1){
  or = order(DNA_CODE$SecondPosition)
  codons = c(rep(DNA_CODE$GeneticCode[or][17:32],3),rep(DNA_CODE$GeneticCode[or][49:64],3))
  tab = as.data.frame(table(sit_mut$mut_to))

  mut_signa = data.frame(codons = codons,
                         to = c(rep("C>A",16),rep("C>T",16),rep("C>G",16),rep("T>A",16),rep("T>C",16),rep("T>G",16)),
                         cancer = c(rep(file,96)))

  mut_signa$codons = factor(mut_signa$codons, levels = c(DNA_CODE$GeneticCode[or][17:32],DNA_CODE$GeneticCode[or][49:64]))
  cod_to = paste0(mut_signa$codons,mut_signa$to)

  mut_signa$counts = countMatches(cod_to,mut_to$mut_to)
  mut_signa$pourcentage = round((mut_signa$counts/sum(mut_signa$counts))*100,2)
  lists = as.character(mut_signa$codons)
  substr(lists,2,2) = "_"
  mut_signa$codons = lists
  return(mut_signa)
}
