#' Count gene set variant
#'
#' @param Table A tibble indicating if each individual has mutation in each of the genes
#' @param GeneSet A one-col table containing the genes to be observed
#' @param Name A string, the name of the gene set, e.g. "high_expressed_genes"
#'
#' @return count_gene_set_variant() returns a table specifying the number of mutations in the genes observed for each and every individual
#' @export
#'
#' @examples 
#' count_gene_set_variant(Table, Set, "high_expressed_genes")

count_gene_set_variant <- function(Table, GeneSet, Name) {
  
  if (is.tibble(Table) == FALSE){
    stop("Table needs to be a tibble.")
  } 
  
  variant.per.gene_set.per.sample<- Table %>%
    dplyr::filter(gene %in% GeneSet) %>%
    dplyr::summarise(across(2:ncol(Table),sum)) %>%
    tidyr::pivot_longer( everything(), names_to = "sampleID", values_to = Name)
  
  return(variant.per.gene_set.per.sample)
}