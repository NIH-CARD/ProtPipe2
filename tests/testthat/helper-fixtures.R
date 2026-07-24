package_path <- function(...) {
  file.path("..", "..", ...)
}

basic_example_path <- function(...) {
  package_path("EXAMPLES", "basic_example_data", ...)
}

load_basic_data <- function() {
  data.table::fread(basic_example_path("iPSC.csv"))
}

load_basic_se <- function() {
  ProtPipe2::create_se(load_basic_data())
}

load_basic_imputed_se <- function() {
  ProtPipe2::impute_min(load_basic_se(), 0)
}
