# Test fixtures are sourced from data the package already ships, so the suite
# runs identically from the source tree and from an installed package (i.e.
# under R CMD check). `neuron_differentiation_intensities` is built verbatim
# from EXAMPLES/basic_example_data/iPSC.csv by
# data-raw/neuron_differentiation_example_data.R, so this is the same data the
# tests previously read off disk.

extdata_path <- function(...) {
  path <- system.file("extdata", ..., package = "ProtPipe2")
  if (!nzchar(path)) {
    testthat::skip(paste0("extdata fixture not installed: ", paste(..., sep = "/")))
  }
  path
}

load_basic_data <- function() {
  data.table::as.data.table(
    get_package_data("neuron_differentiation_intensities")
  )
}

load_basic_se <- function() {
  ProtPipe2::create_se(load_basic_data())
}

load_basic_imputed_se <- function() {
  ProtPipe2::impute_min(load_basic_se(), 0)
}

# utils::data() assigns into an environment rather than returning the object.
get_package_data <- function(name) {
  env <- new.env(parent = emptyenv())
  utils::data(list = name, package = "ProtPipe2", envir = env)
  get(name, envir = env)
}
