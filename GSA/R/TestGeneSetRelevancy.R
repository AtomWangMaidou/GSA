#' Test gene set relevancy
#'
#' @param rawTable A PED format table containing the phenotypes and covariates of each sample.
#' @param GeneSet A one column table containing the col names in the PED table for each gene in the gene set.
#' @param CVnames A one column table containing the col names in the PED table for each covariate in the gene set.
#' @param Trait A character specifying the phenotype to be investigated.
#' @param conf.int -
#'
#' @return A tidy table including regular summary from the glm, e.g. beta, std.error, p.val, odds.ratio, ...
#' @export
#'
#' @examples
#' test_gene_set_relevancy(file.path(R.home("raw_table")), file.path(R.home("gene_set")), file.path(R.home("covariates")), DISEASE)

test_gene_set_relevancy <- function(rawTable, GeneSet, CVnames, Trait, conf.int=FALSE) {
  
  rawTable <- rawTable %>% tibble::as_tibble()
  GeneSet <- GeneSet %>% tibble::as_tibble() %>% pull(1)
  CVnames <- CVnames %>% tibble::as_tibble() %>% pull(1)
  
  glmFunction <- function(inputTable, y){
    f <- base::paste0(y,"~ .")
    stats::glm(formula = as.formula(f), data = inputTable, family = "binomial")
  }
  
  tidy_fun <- function(model, GS_Count, ...) {
    broom::tidy(model, ...) %>%
      dplyr::filter(term == GS_Count)
  }
  
  relevantTable <- rawTable %>% 
    dplyr::select(all_of(GeneSet), all_of(CVnames), all_of(Trait)) %>%
    tidyr::pivot_longer(cols = all_of(GeneSet), names_to = "category", values_to = "count") %>%
    dplyr::group_by(category) %>%
    tidyr::nest() %>%
    dplyr::mutate(logistic_regression = map(data, ~glmFunction(., Trait))) 
  
  relevantTable <- relevantTable %>%
    dplyr::transmute(category, beta = map(logistic_regression, ~tidy_fun(., "count", conf.int = conf.int))) %>%
    tidyr::unnest(cols = c(beta))
  
  return (relevantTable)
  
}