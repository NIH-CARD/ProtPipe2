library(testthat)
setwd("../..")

#spectronaut input
test_that("correctly make a prot_data object", {
  # Print the current working directory
  print(getwd())
  dat <- data.table::fread("EXAMPLES/basic_example_data/iPSC.csv")
  dat1 <- dat[,c(1:2,6:10)]
  dat2 <- dat[,c(1:5)]
  dat_pro_1 <- create_protdata(dat1)
  dat_pro_2 <- create_protdata(dat2)
  dat_pro <- ProtPipe::combine(c(dat_pro_1, dat_pro_2))
  expect_s4_class(dat_pro, "ProtData")
  expect_equal(nrow(dat_pro@condition), 8)
})
