#' Create mutational signature from a cosmic dataset
#'
#' This function will create mutational signature from an input dataset. See each function description
#' for more details.
#'
#' @param file Path to the input file.
#' @param num the number of patient you want to analyse. Default will go through all the dataset.
#' @return change_nu = the change of nucleotide due to the mutation, mut_to = the codon paste to
#' the nucleotide after mutation.
#' @export
Signature_algoV2 = function(file,num, random = FALSE){
  for(i in 1:length(file)){
    if( i == 1){
      print("Running into file 1")
      data = import_data_cos(file[i])
      seq = DNA_sequence_cos(data,num,random)
      mut_to = Site_mut_cos(seq)
      mut_signa = mut_signature_cos(mut_to, i)
      next
    }
    print(c("Running into file", i))
    data = import_data_cos(file[i])
    seq = DNA_sequence_cos(data,num,random)
    mut_to = Site_mut_cos(seq)
    mut_signa1 = mut_signature_cos(mut_to, i)
    mut_signa = rbind(mut_signa,mut_signa1)
  }
  mut_signa$cancer = as.factor(mut_signa$cancer)
  return(mut_signa)
}
