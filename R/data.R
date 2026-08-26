#' Example Proteomics Dataset
#'
#' A small subset of a proteomics dataset used to demonstrate replicate-based
#' imputation.
#'
#' @format A data frame containing protein abundances and experimental
#' conditions in columns. Protein names are set as row names.
#'
#' @source Data derived from a subset of the dataset reported in
#' \doi{10.1016/j.celrep.2025.115362}.
#'
#' @examples
#' data(example_dataset)
#' head(example_dataset)
"example_dataset"


#' Example Proteomics SummarizedExperiment
#'
#' A SummarizedExperiment object containing the example proteomics dataset and
#' experimental design information.
#'
#' @format A \code{SummarizedExperiment} object containing protein abundance
#' measurements and experimental metadata.
#'
#' @source Data derived from a subset of the dataset reported in
#' \doi{10.1016/j.celrep.2025.115362}.
#'
#' @examples
#' data(example_se)
#' example_se
"example_se"
