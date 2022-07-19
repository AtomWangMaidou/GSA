#' Count gene mutations
#'
#' @param table A table containing the information of each and every gene the samples have been tested on.
#' @param gene.col The number of the gene column when counting from the left.
#' @param gt.col The number of the first sample column when counting from the left.
#' @param flag "Allele", "Snps", or "Carrier", indicating the table to be generated
#'
#' @return 
#' flag = Alleles, A table containing the number of mutated alleles for each gene in each sample.
#' flag = Snps, A table containing the number of mutated locus in each gene for each sample.
#' flag = Carrier, A table indicating the presence of mutation in each gene for each sample.
#' 
#' @export
#'
#' @examples
#' count_gene_mutations <- (file.path(R.home("data_table")), 3, 7, Alleles)
#' 
#' count_gene_mutations <- (file.path(R.home("data_table")), 3, 7, Snps)
#' 
#' count_gene_mutations <- (file.path(R.home("data_table")), 3, 7, Carrier)

count_gene_mutations <- function (table, gene.col, gt.col, flag = c("Alleles", "Snps", "Carrier")) {
  
  if(!(is.numeric(gene.col) && gene.col == round(gene.col))){
    stop("gene.col needs to be an integer.")
  }
  
  if(!(is.numeric(gt.col) && gt.col == round(gt.col))){
    stop("gene.col needs to be an integer.")
  } 
  
  flag <- match.arg(flag)
  
  newTable <- dplyr::select(table, c(gt.col:ncol(table)))
  
  gt2count <- function(genotype){
    stringr::str_extract_all(genotype, "\\d") %>%
      purrr::map_int(~sum(as.integer(.)))
  }
  
  topTable <- newTable %>% purrr::map_df(gt2count)
  
  BottomTable <- topTable
  BottomTable[BottomTable >=1 ] <- 1
  
  GeneCol <- dplyr::select(table, c(gene.col))
  
  AllelesRaw <- base::cbind(GeneCol, topTable)
  base::names(AllelesRaw)[1] <- 'gene'
  
  Alleles <- stats::aggregate(. ~ gene, data=AllelesRaw, FUN=sum)
  
  SnpsCarrierRaw <- base::cbind(GeneCol, BottomTable)
  base::names(SnpsCarrierRaw)[1] <- 'gene'
  
  Snps <- stats::aggregate(. ~ gene, data=SnpsCarrierRaw, FUN=sum)
  
  Carrier <- Snps
  Carrier[Carrier >=1 ] <- 1
  
  if (flag == "Alleles"){
    return (Alleles)
  }
  else if (flag == "Snps"){
    return (Snps)}
  else if (flag == "Carrier"){
    return (Carrier)
  }
  
}




