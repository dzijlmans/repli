#' Reconstruct a SummarizedExperiment after replicate-based imputation
#'
#' Places imputed data from [impute_replicate()] back into the original
#' SummarizedExperiment object supplied in [prepare_imputation_data_se()].
#'
#' @param imputed_data A data frame containing the imputed data returned
#'   by [impute_replicate()].
#'
#' @param se The original SummarizedExperiment object supplied in
#' [prepare_imputation_data_se()].
#'
#' @param experimental_design A data frame describing the experimental
#'   design, containing columns \code{label}, \code{condition}, and
#'   \code{replicate}.
#'
#' @returns A SummarizedExperiment object containing the imputed assay. The
#' original row and column metadata are retained.
#'
#' @export
#'
#' @examples
#' data(example_se)
#' prepared <- prepare_imputation_data_se(example_se)
#' imputed <- impute_replicate(prepared)
#' result <- reconstruct_se(
#'   imputed,
#'   example_se,
#'   prepared$experimental_design
#' )
reconstruct_se <- function(imputed_data, se, experimental_design) {

  df <- imputed_data %>%
    dplyr::select(-sd, -row_mean, -missingness) %>%
    tidyr::pivot_longer(
      cols = -c(name, condition),
      names_to = "replicate",
      values_to = "value"
    ) %>%
    dplyr::left_join(experimental_design, by = c("condition", "replicate")) %>%
    dplyr::select(name, label, value) %>%
    tidyr::pivot_wider(names_from = label, values_from = value)

  df <- df %>%
    dplyr::slice(match(rownames(se), name))

  value_cols <- experimental_design$label

  imputed_matrix <- as.matrix(df[, value_cols, drop = FALSE])

  rownames(imputed_matrix) <- df$name
  colnames(imputed_matrix) <- value_cols

  SummarizedExperiment::assay(se) <- imputed_matrix

  se
}
