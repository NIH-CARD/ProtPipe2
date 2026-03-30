test_that("numeric utility helpers behave on mixed input data", {
  df <- data.frame(
    gene = c("A", "B"),
    chr_num = c("1.5", "2.5"),
    chr_text = c("x", "y"),
    num = c(10, 20),
    stringsAsFactors = FALSE
  )

  converted <- ProtPipe::convert_numeric_cols(df)

  expect_true(is.numeric(converted$chr_num))
  expect_true(is.character(converted$chr_text))
  expect_equal(unname(ProtPipe::detect_intensity_cols(converted)), c(2, 4))
})

test_that("ProtData sample counting works on bundled example data", {
  dat <- load_basic_data()
  pd <- ProtPipe::create_protdata(dat)

  expect_s4_class(pd, "ProtData")
  expect_equal(ProtPipe::num_samples(pd), 42)
})
