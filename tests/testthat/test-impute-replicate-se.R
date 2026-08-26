test_that("impute_replicate_se preserves SummarizedExperiment structure",
          {
            se <- test_se()

            set.seed(123)
            result <- impute_replicate_se(se, p = 1)

            expect_s4_class(result, "SummarizedExperiment")

            expect_identical(dim(result), dim(se))

            expect_identical(rownames(result), rownames(se))

            expect_identical(colnames(result), colnames(se))
          })

test_that("impute_replicate_se imputes MAR values", {
  se <- test_se()

  set.seed(123)
  result <- impute_replicate_se(se, p = 1)

  original <- SummarizedExperiment::assay(se)
  imputed <- SummarizedExperiment::assay(result)

  expect_true(!is.na(imputed["protein_2", "control_2"]))

  expect_true(!is.na(imputed["protein_3", "treatment_2"]))

  expect_true(all(is.na(imputed["protein_4", 4:6])))
})
