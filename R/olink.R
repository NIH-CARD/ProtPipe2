#' Create a ProtData Object from Olink NPX Data
#'
#' @description
#' A wrapper function that processes a data frame of Olink NPX data, optionally
#' merges it with a sample condition/metadata file, and formats the result into a
#' ProtData object for downstream analysis.
#'
#' @param npx A data frame containing Olink NPX data, typically from OlinkAnalyze::read_npx().
#' @param condition An optional data frame containing sample metadata. It is crucial
#'   that this file contains a column named 'SampleID' for merging.
#' @param filter A logical value (TRUE/FALSE) indicating whether to filter out
#'   samples with QC warnings. Defaults to FALSE.
#'
#' @return A ProtData object formatted for use with other package functions.
#'
#' @export
#'
#' @importFrom dplyr select everything
#'
create_protdata_from_olink <- function(npx, condition = NULL, filter = F) {
  npx <- as.data.frame(npx)
  dat <- olink_sample_out(npx, filter)
  if(!is.null(condition)){
    if (!"SampleID" %in% colnames(condition)) {
      stop("Error: The uploaded file must contain a column named 'SampleID'.")
    }
    condition <- dplyr::select(condition, SampleID, dplyr::everything())
  }
  return(create_protdata(dat, intensity_cols = c(4:length(colnames(dat))), condition, method = "Olink"))
}


olink_sample_out=function(my_npx, filter){
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
