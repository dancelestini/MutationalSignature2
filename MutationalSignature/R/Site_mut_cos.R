#' Substitute same molecular mutation
#'
#' This function will just change nucleotide from mutation to get only T and C mutated.
#'
#' @param seq Output from DNA_sequence_cos function
#' @return change_nu = the change of nucleotide due to the mutation, mut_to = the codon paste to
#' the nucleotide after mutation.
#' @export
Site_mut_cos = function(seq){
  seq2 = seq
  mut_to = c()
  for( i in 1:length(seq2$change_nu)){
    from = strsplit(seq2$change_nu[i], "")[[1]][1]
    to = strsplit(seq2$change_nu[i], "")[[1]][3]
    if(from == "A"){# To invert all change involving Adenosine
      from = "T"
      if (to == "T"){
        to = "A"
      }
      if (to == "C"){
        to = "G"
      }
      if (to == "G"){
        to = "C"
      }
    }
    if(from == "G"){# To invert all change involving Guanine
      from = "C"
      if (to == "T"){
        to = "A"
      }
      if (to == "C"){
        to = "G"
      }
      if (to == "A"){
        to = "T"
      }
    }
    if (!is.na(seq2$codon_part[[i]][3])){
      seq2$codon_part[[i]][2] = from
      seq2$change_nu[i] = paste0(from,">",to)
      mut_to = c(mut_to,paste0(seq2$codon_part[[i]][1],seq2$codon_part[[i]][2],seq2$codon_part[[i]][3],seq2$change_nu[i]))
      next
    }
    mut_to = c(mut_to, NA_character_)
    seq2$change_nu[i] = NA_character_
  }
  return(list(change_nu = seq2$change_nu,mut_to = mut_to))
}
