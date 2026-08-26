test_that("mixed imputation leaves MNAR values missing when MNAR method is none",
          {
            skip_if_not_installed("DEP")

            se <- test_se()

            set.seed(123)

            result <- impute_replicate_mixed_se(se, p = 1, MNAR_method = "none")

            assay_result <- SummarizedExperiment::assay(result)

            # MAR is imputed
            expect_true(!is.na(assay_result["protein_2", "control_2"]))

            expect_true(!is.na(assay_result["protein_3", "treatment_2"]))

            # MNAR remains missing
            expect_true(all(is.na(assay_result["protein_4", 4:6])))
          })
