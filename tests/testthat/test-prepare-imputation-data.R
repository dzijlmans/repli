test_that("prepare_imputation_data prepares data correctly", {
  data <- test_data()

  prepared <- prepare_imputation_data(data, p = 1)

  expect_named(prepared, c("data", "replicate_cols", "global_sd"))

  expect_equal(length(prepared$replicate_cols), 3)

  expect_true(is.numeric(prepared$global_sd))
})


test_that("prepare_imputation_data correctly classifies missingness", {
  prepared <- prepare_imputation_data(test_data())

  result <- prepared$data

  expect_true(all(result$missingness[result$name == "protein_1"] == "Valid"))

  expect_true(any(result$missingness[result$name == "protein_2" &
                                       result$condition == "control"] == "MAR"))

  expect_true(any(result$missingness[result$name == "protein_3" &
                                       result$condition == "treatment"] == "MAR"))

  expect_true(any(result$missingness[result$name == "protein_4" &
                                       result$condition == "treatment"] == "MNAR"))
})
