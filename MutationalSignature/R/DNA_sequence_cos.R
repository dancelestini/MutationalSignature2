 #' Obtain all mutations codons
#'
#' This function will load all the CDS from the mutated genes in the input file. Then recovery the 3' and
#' 5' nucleotide aside of the mutation site.
#'
#' @param cosmic Dataframe output from the import_data_cos function
#' @param num The number of patient you want to analyse. Default will go through all the dataset.
#' @return DNA_seq = CDS sequence for all the genes used, codon_part = codon around mutation.
#' @export
DNA_sequence_cos = function(cosmic,num = length(cosmic$GENE_NAME),random = FALSE){
  if (random == TRUE){
    idx = sort(sample(1:length(cosmic$GENE_NAME),num, replace = F))
    cosmic = cosmic[idx,]
  }
  genes_mut_names = sapply(strsplit(cosmic$GENE_NAME[1:num],"[\\_]+"), `[`, 1)
  genes_mut_site1 = sapply(strsplit(cosmic$MUTATION_CDS[1:num],"[\\c.]+"), `[`, 2)
  genes_mut_site = as.numeric(str_extract(genes_mut_site1, "[0-9]+"))
  genes_mut_site2 = sapply(strsplit(cosmic$MUTATION_AA[1:num],"[\\p.]+"), `[`, 2)
  change_nu = sapply(strsplit(genes_mut_site1,"[0-9]+"),tail,1)
  change_aa =substr(genes_mut_site2, nchar(genes_mut_site2), nchar(genes_mut_site2))
  change_aa[change_aa %ni% DNA_CODE$AA] = NA_character_

  term_genes = paste0(rep("Human[orgn] AND ",length(genes_mut_names)),genes_mut_names, rep("[gene]",length(genes_mut_names)))
  same_genes = as.numeric(factor(term_genes, levels = unique(term_genes)))

  DNA_seq = c()
  codon_part = c()
  name_seq = c()
  print(paste0("Number of genes to load = ",length(unique(term_genes))))
  for (i in 1:length(unique(term_genes))){
    print(i)
    uids =esearch(db = "gene",term = unique(term_genes)[i],retmax = 5)
    x = elink(uids, dbFrom="gene", dbTo="nuccore")
    x =linkset(x)
    # if(is.na(uids[1]) || is.na(change_aa[i])){
    #   DNA_seq = c(DNA_seq,list(NA_character_))
    #   codon_part = c(codon_part,list(NA_character_))
    #   next
    # }
    access = efetch(x$gene_nuccore[1:5], rettype = "acc",retmode = "text")
    access = strsplit(content(access), "\n")[[1]]
    access = access[grep("NM_",access)]
    if(is.na(access[1])){
      DNA_seq = c(DNA_seq,list(NA_character_))
      next
    }
    longest = FindLongestSeq(access)$Accession
    access = efetch(db = "nuccore",uid = longest, rettype = "fasta_cds_na", retmode = "text")
    access = strsplit(content(access), "\n")[[1]]

    if(access[1] ==""){
      DNA_seq = c(DNA_seq,list(NA_character_))
      next
    }
    name_seq = c(name_seq,access[1])
    access = access[-1]
    access =sapply(strsplit(access,">"), `[`, 1)
    idx = which(access == "" | is.na(access))
    access = access[1:idx[1]-1]
    DNA_seq_part = c()
    for (j in access){
      DNA_seq_part = c(DNA_seq_part,strsplit(j,"")[[1]])
    }
    DNA_seq = c(DNA_seq,list(DNA_seq_part))
  }

  for(k in 1:length(term_genes)){
    if(is.na(DNA_seq[[same_genes[k]]][1]) || is.na(genes_mut_site[k]) || is.na(DNA_seq[[same_genes[k]]][(genes_mut_site[k]-1):(genes_mut_site[k]+1)])){
      codon_part = c(codon_part,list(NA_character_))
      print(k)
      next
    }

    codon_part = c(codon_part,list(DNA_seq[[same_genes[k]]][(genes_mut_site[k]-1):(genes_mut_site[k]+1)]))
  }

  return(list(DNA_seq = DNA_seq, codon_part = codon_part,genes_mut_names = genes_mut_names,
              genes_mut_site = genes_mut_site, change_aa = change_aa,change_nu = change_nu))
}
