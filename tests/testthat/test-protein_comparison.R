dat <- data.table::fread("EXAMPLES/basic_example_data/iPSC.csv")
dat_pro <- create_se(dat)

test_that("correctly compare proteins", {
  p <- ProtPipe::compare_protein(dat_pro, "U3KQP1")
  p <- ProtPipe::compare_protein(dat_pro, "U3KQP1", condition = "base_condition")
  p <- ProtPipe::compare_protein(dat_pro, "U3KQP1", condition = "base_condition",
                                 selected_groups = c("Day0", "Day28") )
})
