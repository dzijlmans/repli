test_that("impute_replicate imputes MAR values and leaves observed values unchanged", {
  prepared <- prepare_imputation_data(test_data(), p = 1)

  set.seed(123)
  result <- impute_replicate(prepared)

  original <- as.matrix(prepared$data[, prepared$replicate_cols])

  imputed <- as.matrix(result[, prepared$replicate_cols])

  mar <- prepared$data$missingness == "MAR"
  missing <- is.na(original)

  expect_true(all(!is.na(imputed[mar & missing])))

  observed <- !is.na(original)

  expect_equal(imputed[observed], original[observed])
})

test_that("impute_replicate does not impute MNAR values", {
  prepared <- prepare_imputation_data(test_data(), p = 1)

  set.seed(123)
  result <- impute_replicate(prepared)

  mnar <- prepared$data$missingness == "MNAR"

  expect_true(all(is.na(as.matrix(result[mnar, prepared$replicate_cols]))))
})
