package_path <- function(...) {
  file.path("..", "..", ...)
}

basic_example_path <- function(...) {
  package_path("EXAMPLES", "basic_example_data", ...)
}

load_basic_data <- function() {
  data.table::fread(basic_example_path("iPSC.csv"))
}

load_basic_metadata <- function(include_numeric = FALSE, include_batch = FALSE) {
  se <- ProtPipe::create_se(load_basic_data())
  meta <- as.data.frame(SummarizedExperiment::colData(se))
  meta$SampleID <- rownames(meta)

  if (include_numeric) {
    meta$day_num <- as.numeric(gsub("\\D", "", meta$base_condition))
  }

  if (include_batch) {
    meta$batch <- rep(c("Batch1", "Batch2"), length.out = nrow(meta))
  }

  meta[, c("SampleID", setdiff(names(meta), "SampleID")), drop = FALSE]
}

load_basic_se <- function(sample_metadata = NULL) {
  ProtPipe::create_se(load_basic_data(), sample_metadata = sample_metadata)
}

load_basic_imputed_se <- function(sample_metadata = NULL) {
  ProtPipe::impute_min(load_basic_se(sample_metadata), 0)
}
