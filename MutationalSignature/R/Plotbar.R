#' Create mutational signature from a cosmic dataset
#'
#' This function will create mutational signature from an input dataset. See each function description
#' for more details.
#'
#' @param total Path to the output of Signature_algoV2 function.
#' @return Mutational signature using ggplot2.
#' @export
Plotbar = function(total){
  ggplot(total,aes(x=codons,y=pourcentage, fill = to))+
    geom_bar(stat = "identity", size = 0.5, width = 0.8)+
    ggtitle("Cosmic dataset mutational signature")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(size=14, face="bold.italic")) +
    facet_grid(cols = vars(to))
}
