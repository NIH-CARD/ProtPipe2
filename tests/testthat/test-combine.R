#spectronaut input
test_that("correctly make a prot_data object", {
  dat <- data.table::fread(basic_example_path("iPSC.csv"))
  dat1 <- dat[,c(1:2,6:10)]
  dat2 <- dat[,c(1:5)]
  dat_pro_1 <- create_protdata(dat1)
  dat_pro_2 <- create_protdata(dat2)
  dat_pro <- ProtPipe::combine(c(dat_pro_1, dat_pro_2))
  expect_s4_class(dat_pro, "ProtData")
  expect_equal(nrow(dat_pro@condition), 8)
})
