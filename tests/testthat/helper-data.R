test_data <- function() {
  data <- matrix(
    c(10, 10.5, 11, 12, 12.5, 13,
      10, NA, 11, 12, 12.5, 13,
      10, 10.5, 11, 12, NA, 13,
      10, 10.5, 11, NA, NA, NA,
      20, 20.5, 21, 18, 18.5, 19),
    nrow = 5,
    byrow = TRUE
  )


  rownames(data) <- paste0("protein_", 1:5)

  colnames(data) <- c("control_1",
                      "control_2",
                      "control_3",
                      "treatment_1",
                      "treatment_2",
                      "treatment_3")

  as.data.frame(data)
}

test_se <- function() {
  data <- test_data()

  data <- 2^data

  sim <- data %>%
    tibble::rownames_to_column("name") %>%
    dplyr::mutate(ID = seq_len(dplyr::n()))

  experimental_design <- data.frame(
    label = colnames(data),
    condition = c(
      "control",
      "control",
      "control",
      "treatment",
      "treatment",
      "treatment"
    ),
    replicate = c(1, 2, 3, 1, 2, 3)
  )

  value_cols <- which(
    colnames(sim) %in% experimental_design$label
  )

  DEP::make_se(sim, value_cols, experimental_design)
}
