#' Create a SummarizedExperiment Object from Olink NPX Data
#'
#' @description
#' A wrapper function that processes a data frame of Olink NPX data, optionally
#' merges it with a sample condition/metadata file, and formats the result into a
#' SummarizedExperiment object for downstream analysis.
#'
#' @param npx A data frame containing Olink NPX data, typically from OlinkAnalyze::read_npx().
#' @param condition An optional data frame containing sample metadata. It is crucial
#'   that this file contains a column named 'SampleID' for merging.
#' @param filter A logical value (TRUE/FALSE) indicating whether to filter out
#'   samples with QC warnings. Defaults to FALSE.
#'
#' @return A SummarizedExperiment object formatted for use with other package functions.
#'
#' @export
#'
#' @importFrom dplyr select everything
#'
create_se_from_olink <- function(npx, condition = NULL, filter = TRUE) {
  protpipe_require_packages("tidyr", feature = "create_se_from_olink()")

  npx <- as.data.frame(npx)
  dat <- olink_sample_out(npx, filter)
  if(!is.null(condition)){
    if (!"SampleID" %in% colnames(condition)) {
      stop("Error: The uploaded file must contain a column named 'SampleID'.")
    }
    condition <- dplyr::select(condition, SampleID, dplyr::everything())
  }
  # create_se() marks Olink objects as already log2 scaled via creation_method.
  return(create_se(data = dat, intensity_cols = c(4:length(colnames(dat))),
                   sample_metadata = condition, creation_method = "Olink"))
}


#' Format Olink NPX Data into Data and Condition Tables
#'
#' @param my_npx A data frame of Olink NPX measurements.
#'
#' @return A list containing `data`, `condition`, and `number_samples`.
#' @export
olink_all_output <- function(my_npx){
  protpipe_require_packages("tidyr", feature = "olink_all_output()")

  # Compute per-protein max LOD across plates before pivoting to avoid duplicate rows
  lod_per_protein <- my_npx |>
    dplyr::group_by(UniProt, Assay, OlinkID) |>
    dplyr::summarise(LOD = max(LOD, na.rm = TRUE), .groups = "drop")

  npx_wide <- my_npx |>
    dplyr::select(SampleID, UniProt, Assay, OlinkID, NPX) |>
    tidyr::pivot_wider(names_from = SampleID, values_from = NPX, values_fn = mean) |>
    dplyr::rename(Protein_Group = UniProt, Genes = Assay) |>
    dplyr::left_join(
      lod_per_protein |> dplyr::rename(Protein_Group = UniProt, Genes = Assay),
      by = c("Protein_Group", "Genes", "OlinkID")
    ) |>
    dplyr::select(Protein_Group, Genes, OlinkID, LOD, dplyr::everything()) |>
    as.data.frame()

  number_samples <- ncol(npx_wide) - 4
  condition <- NULL
  return(list(data = npx_wide, condition = condition, number_samples = number_samples))
}

olink_sample_out=function(my_npx, filter = T){
  protpipe_require_packages("tidyr", feature = "olink_sample_out()")

  if(filter){
    npx_wide <- my_npx |>
      dplyr::mutate(NPX = ifelse(NPX < LOD, NA, NPX)) |>
      # dplyr::filter(SampleType == "SAMPLE") |>
      # dplyr::filter(AssayType == "assay") |>
      dplyr::select(SampleID, UniProt,Assay, OlinkID, NPX) |>
      tidyr::pivot_wider(names_from = SampleID, values_from = NPX, values_fn = mean) |>
      dplyr::rename(Protein_Group = UniProt, Genes = Assay)
  }else{
    npx_wide <- my_npx |>
      # dplyr::filter(AssayType == "assay") |>
      dplyr::select(SampleID, UniProt,Assay, OlinkID, NPX) |>
      tidyr::pivot_wider(names_from = SampleID, values_from = NPX,values_fn = mean) |>
      dplyr::rename(Protein_Group = UniProt, Genes = Assay)
  }
  return(npx_wide)
}
