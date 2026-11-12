#' Perform mixed MAR/MNAR imputation on a SummarizedExperiment object
#'
#' Performs mixed imputation of missing values in a SummarizedExperiment object.
#' Missing At Random (MAR) values are imputed using replicate-based imputation,
#' while Missing Not At Random (MNAR) values are imputed using a specified
#' method from the DEP R package. Both methods are applied independently to
#' the original dataset, after which their results are combined.
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
#' @param MNAR_method A character string specifying the DEP imputation method
#' used for MNAR values, passed to [DEP::impute()]. Default is \code{"MinProb"}.
#'
#' @param MNAR_args A named list of additional arguments passed to
#' [DEP::impute()]. Default is \code{list(q = 0.01)}.
#'
#' @details Protein-condition pairs are classified as MAR if there are missing
#' values for some but not all replicates. Protein-condition pairs are
#' classified as MNAR if there are missing values for all of the replicates of
#' a condition.
#' 
#' #' MAR and MNAR imputation are performed independently on the original dataset.
#' Missing MAR values are imputed using replicate-based imputation, while
#' missing MNAR values are imputed using the specified DEP method. The results
#' are subsequently combined, retaining all originally observed values.
#'
#' When \code{MNAR_method = "none"}, MNAR values are left missing.
#'
#' @returns A SummarizedExperiment object containing the imputed assay. MAR
#' values are imputed using replicate-based imputation and MNAR values are
#' imputed using the specified DEP method.
#'
#' @export
#'
#' @examples
#' data(example_se)
#' imputed_se <- impute_replicate_mixed_se(
#'   example_se,
#'   p = 0.05,
#'   MNAR_method = "MinProb",
#'   MNAR_args = list(q = 0.01)
#' )
impute_replicate_mixed_se <- function(se,
                                      p = 0.05,
                                      MNAR_method = "MinProb",
                                      MNAR_args = list(q = 0.01)) {

  if (!requireNamespace("DEP", quietly = TRUE)) {
    stop(
      "The 'DEP' package is required for impute_replicate_mixed_se(). ",
      "Please install it first."
    )
  }

  df_long <- DEP::get_df_long(se)

  proteins_MNAR <- df_long %>%
    dplyr::group_by(name, condition) %>%
    dplyr::summarise(
      all_missing = all(is.na(intensity)),
      .groups = "drop"
    ) %>%
    dplyr::filter(all_missing) %>%
    dplyr::pull(name) %>%
    unique()

  X <- SummarizedExperiment::assay(se)
  MNAR_rows <- rownames(X) %in% proteins_MNAR
  mnar_matrix <- is.na(X) & MNAR_rows
  mar_matrix <- is.na(X) & !mnar_matrix

  se_MAR <- impute_replicate_se(
    se,
    p = p
  )
  
  if (MNAR_method != "none") {
    
    se_MNAR <- do.call(
      DEP::impute,
      c(list(se, fun = MNAR_method), MNAR_args)
    )
    
  }
  
  X_final <- X
  
  X_MAR <- SummarizedExperiment::assay(se_MAR)
  X_final[mar_matrix] <- X_MAR[mar_matrix]
  
  if (MNAR_method != "none") {
    
    X_MNAR <- SummarizedExperiment::assay(se_MNAR)
    X_final[mnar_matrix] <- X_MNAR[mnar_matrix]
    
  }
  
  se_final <- se
  SummarizedExperiment::assay(se_final) <- X_final
  
  return(se_final)
}
