#' Perform replicate-based imputation
#'
#' Performs replicate-based imputation of Missing At Random (MAR) values.
#' Missing values are imputed by sampling from a normal distribution using the
#' mean of the observed replicates for the corresponding protein-condition pair
#' and a global standard deviation derived from replicates of all
#' protein-condition pairs across the dataset.
#'
#' @param prepared_data An object returned by
#' [prepare_imputation_data()] or [prepare_imputation_data_se()]. Contains a
#' list consisting of data in a format prepared for imputation (including
#' missingness classification and protein-condition means), a vector specifying
#' the columns containing replicate measurements and the global standard
#' deviation.
#'
#' @returns A data frame containing imputed data. MAR values are imputed,
#' while values classified as MNAR or observed are left unchanged.
#'
#' @export
#'
#' @examples
#' #in this instance data is a data frame with proteins as rows and
#' #experimental conditions as columns. Protein names are set as rownames.
#' data(example_dataset)
#' prepared <- prepare_imputation_data(example_dataset)
#' imputed <- impute_replicate(prepared)
impute_replicate <- function(prepared_data) {

  data_long <- prepared_data$data
  replicate_cols <- prepared_data$replicate_cols
  global_sd <- prepared_data$global_sd

  mat <- as.matrix(data_long[, replicate_cols])

  mar_mask <- data_long$missingness == "MAR"
  na_mask <- is.na(mat)

  for (j in seq_len(ncol(mat))) {
    rows_to_impute <- which(mar_mask & na_mask[, j])

    if (length(rows_to_impute) > 0) {
      mat[rows_to_impute, j] <- stats::rnorm(length(rows_to_impute),
                                             mean = data_long$row_mean[rows_to_impute],
                                             sd = global_sd)
    }
  }

  data_long[, replicate_cols] <- as.data.frame(mat)

  data_long
}
