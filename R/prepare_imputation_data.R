#' Prepare data for replicate-based imputation
#'
#' Prepares a data frame containing protein abundance values for
#' replicate-based imputation. Data is transformed in long format, with
#' protein-condition pairs as rows and replicates as columns. Missingness is
#' classified as Missing At Random (MAR), Missing Not At Random (MNAR) or
#' observed. For protein-condition pairs classified as MAR, the mean abundance
#' is calculated. A global standard deviation is derived from the standard
#' deviations of replicates of all protein-condition pairs across the dataset.
#'
#' @param data A data frame containing protein abundance values, with protein
#' names as row names and experimental conditions as column names. Column names
#' should follow condition-replicate naming convention such as
#' \code{"control_1"} or \code{"treatment.2"}. An underscore or period must
#' separate the condition and replicate identifier.
#'
#' @param p A numeric value between 0 and 1 specifying the proportion of
#' proteins with the highest standard deviation across replicates, used to
#' estimate the global standard deviation. Default is \code{0.05}. This means
#' the global standard deviation will be calculated as the median of the top 5%
#' highest standard deviations across all protein-condition pairs in the
#' dataset.
#'
#' @details Protein-condition pairs are classified as MAR if there are missing
#' values for some but not all replicates. Protein-condition pairs are
#' classified as MNAR if there are missing values for all of the replicates of
#' a condition.
#'
#' @returns A list containing:
#' \describe{
#'   \item{data}{The data in a format prepared for imputation, including
#'   missingness classification and protein-condition means.}
#'   \item{replicate_cols}{A vector specifing the replicate columns.}
#'   \item{global_sd}{The global standard deviation among the proteins
#'   in the upper \code{p} proportion of replicate variability.}
#' }
#'
#' @export
#'
#' @examples
#' data(example_dataset)
#' prepared <- prepare_imputation_data(example_dataset)
prepare_imputation_data <- function(data, p = 0.05) {

  data_long <- data %>%
    tibble::rownames_to_column("name") %>%
    tidyr::pivot_longer(
      cols = -name,
      values_to = "values",
      names_to = c("condition", "replicate"),
      names_pattern = "^(.*)[_.](\\d+)$"
    ) %>%
    tidyr::pivot_wider(names_from = "replicate", values_from = "values")

  replicate_cols <- setdiff(names(data_long), c("name", "condition"))

  mat <- as.matrix(data_long[, replicate_cols])

  data_long$sd <- matrixStats::rowSds(mat, na.rm = TRUE)

  median_top_sd <- data_long %>%
    dplyr::arrange(dplyr::desc(sd)) %>%
    dplyr::slice_head(prop = p) %>%
    dplyr::summarise(median_var = stats::median(sd, na.rm = TRUE)) %>%
    dplyr::pull(median_var)

  all_na <- rowSums(!is.na(mat)) == 0
  any_na <- rowSums(is.na(mat)) > 0

  data_long$missingness <- ifelse(all_na, "MNAR", ifelse(any_na, "MAR", "Valid"))

  row_means <- rowMeans(mat, na.rm = TRUE)
  row_means[data_long$missingness != "MAR"] <- NA

  data_long$row_mean <- row_means

  list(
    data = data_long,
    replicate_cols = replicate_cols,
    global_sd = median_top_sd
  )
}
