test_that("correctly perform ttest", {
  df <- readRDS("EXAMPLES/VIRUS/virus_data.rds")
  meta <- readRDS("EXAMPLES/VIRUS/virus_metadata.rds")

  dat_pro <- ProtPipe::create_protdata(df, condition = meta)
  ProtPipe::plot_CVs(dat_pro, "viral.exposure")
  ProtPipe::plot_CVs_optim(dat_pro, "viral.exposure")
  ProtPipe::plot_CVs_optim(dat_pro, "viral.exposure", plot_type = "jitter")
  system.time({
    ProtPipe::plot_CVs(dat_pro, "viral.exposure")
  })
  system.time({
    ProtPipe::plot_CVs_optim(dat_pro, "viral.exposure")
  })

})
