#' Perform replicate-based imputation on a SummarizedExperiment object
#'
#' Performs replicate-based imputation of Missing At Random (MAR) values
#' in a SummarizedExperiment object. Missing values are imputed by sampling
#' from a normal distribution using the mean of the observed replicates for the
#' corresponding protein-condition pair and a global standard deviation derived
#' from replicates of all protein-condition pairs across the dataset.The assay
#' is prepared using [prepare_imputation_data_se()], imputation is performed
#' using [impute_replicate()], and the imputed values are returned in the
#' original SummarizedExperiment structure using [reconstruct_se()].
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
#' @returns A SummarizedExperiment object with MAR values imputed. Completely
#' missing protein-condition pairs, classified as MNAR, are not imputed.
#'
#' @export
#'
#' @examples
#' data(example_se)
#' imputed_se <- impute_replicate_se(example_se)
impute_replicate_se <- function(se, p = 0.05) {

  prepared <- prepare_imputation_data_se(se = se, p = p)

  imputed <- impute_replicate(prepared)

  reconstruct_se(
    imputed_data = imputed,
    se = se,
    experimental_design = prepared$experimental_design
  )
}
