#' Prepare a SummarizedExperiment for replicate-based imputation
#'
#' Extracts abundance data from a SummarizedExperiment and prepares it
#' for replicate-based imputation using [prepare_imputation_data()].
#' The experimental design is also extracted from the sample names for
#' downstream reconstruction of the SummarizedExperiment.
#'
#' @param se A \code{SummarizedExperiment} object containing protein abundance
#' measurements. Sample names must contain a condition and replicate identifier
#' separated by an underscore or period, such as \code{"control_1"} or
#' \code{"treatment.2"}.
#'
#' @param p A numeric value between 0 and 1 specifying the proportion of
#' proteins with the highest standard deviation across replicates, used to
#' estimate the global standard deviation. Default is \code{0.05}. This means
#' the global standard deviation will be calculated as the median of the top 5%
#' highest standard deviations across all protein-condition pairs in the
#' dataset.
#'
#' @returns A list containing:
#' \describe{
#'   \item{data}{The data in a format prepared for imputation, including
#'   missingness classification and protein-condition means.}
#'   \item{replicate_cols}{A vector specifing the replicate columns.}
#'   \item{global_sd}{The global standard deviation among the proteins
#'   in the upper \code{p} proportion of replicate variability.}
#'   \item{experimental_design}{A data frame describing the experimental
#'   design, containing columns \code{label}, \code{condition}, and
#'   \code{replicate}}
#' }
#'
#' @export
#'
#' @examples
#' data(example_se)
#' prepared <- prepare_imputation_data_se(example_se)
#' imputed <- impute_replicate(prepared)
prepare_imputation_data_se <- function(se, p = 0.05) {

  data <- SummarizedExperiment::assay(se) %>%
    as.data.frame()

  prepared_data <- prepare_imputation_data(data = data, p = p)

  prepared_data$experimental_design <- data.frame(label = colnames(data)) %>%
    tidyr::extract(
      label,
      into = c("condition", "replicate"),
      regex = "^(.*)[_.](.*)$",
      remove = FALSE
    )

  prepared_data

}
