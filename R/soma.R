#' Create a SummarizedExperiment Object from SomaScan Data
#'
#' @description
#' Processes SomaScan data from an ADAT object, optionally performs filtering
#' (e.g., on buffer controls), and formats the result into a SummarizedExperiment object
#' for downstream analysis.
#'
#' @param adat A SomaScan ADAT object, typically loaded using the SomaDataIO package.
#' @param condition An optional data frame containing sample metadata. This is passed to
#'   the downstream create_protdata function for integration.
#' @param filter A logical value (TRUE/FALSE) indicating whether to perform
#'   filtering on the data, including for buffer controls. Defaults to TRUE.
#'
#' @return A SummarizedExperiment object formatted for use with other package functions.
#'
#' @export
#'
#' @importFrom dplyr rename
#' @importFrom magrittr %>%
#'
create_se_from_soma <- function(adat, condition = NULL, filter = TRUE) {

  soma_out <- soma_all_output(adat)
  dat <- soma_out$data
  soma_condition <- soma_out$condition
  number_samples <- soma_out$number_samples

  if(filter){
    dat <- Buffer_filter(dat)
  }

  # combine condition with soma sample metadata
  if(!is.null(condition)){
    condition$SampleID <- as.character(condition$SampleID)
    soma_condition <- soma_condition %>%
      dplyr::left_join(condition, by = "SampleID")
  }

  # get intensity cols
  num_rows <- length(names(dat))
  intensity_cols <- c((num_rows-number_samples+1):num_rows)

  return(create_se(data = dat, sample_metadata = soma_condition, intensity_cols = intensity_cols, creation_method = "SomaScan"))
}

# soma_sample_out=function(DT){
#   anno <- SomaDataIO::getAnalyteInfo(DT)%>%
#     dplyr::filter(Organism == "Human") %>%
#     dplyr::filter(Type == "Protein")
#   DT=as.data.frame(DT)
#   DT_dat=DT%>%
#     dplyr::filter(grepl("Sample", SampleType, ignore.case = TRUE))
#
#   #check for duplicated SampleID
#   duplicate_ids <- DT_dat$SampleId[duplicated(DT_dat$SampleId) | duplicated(DT_dat$SampleId, fromLast = TRUE)]
#   if(length(duplicate_ids>0)){
#     cat(paste0("removing duplicates: ", paste(duplicate_ids, collapse = ", ")))
#     DT_dat <- DT_dat[!DT_dat$SampleId %in% duplicate_ids, ]
#   }
#
#
#   rownames(DT_dat)=DT_dat$SampleId
#   condition=DT_dat %>%
#     dplyr::select(-matches("seq\\.", ignore.case = TRUE))
#   DT_dat=DT_dat%>%
#     dplyr::select(matches("seq\\.", ignore.case = TRUE))%>%
#     t()
#   DT_dat=merge(anno[,grep('AptName|UniProt|EntrezGeneSymbol|TargetFullName',colnames(anno))], DT_dat,by.x='AptName',by.y=0,all.x=T)
#   DT_dat=DT_dat %>%
#     dplyr::filter(UniProt != "")%>%
#     dplyr::filter(EntrezGeneSymbol != "") %>%
#     dplyr::rename(Protein_Group= UniProt)%>%
#     dplyr::rename(Genes= EntrezGeneSymbol)
#   return(list(data = DT_dat, condition = condition))
# }


#' format soma adat into data and condition dataframes
#'
#' @param DT 
#'
#' @return
#' @export
#'
#' @examples
soma_all_output <- function(DT){
  anno=SomaDataIO::getAnalyteInfo(DT)
  DT=data.frame(DT)
  DT_dat=data.frame(DT)%>%
    dplyr::filter(grepl("Sample", SampleType, ignore.case = TRUE))


  #check for duplicated SampleID
  duplicate_ids <- DT_dat$SampleId[duplicated(DT_dat$SampleId) | duplicated(DT_dat$SampleId, fromLast = TRUE)]
  if(length(duplicate_ids>0)){
    cat(paste0("removing duplicates: ", paste(duplicate_ids, collapse = ", ")))
    DT_dat <- DT_dat[!DT_dat$SampleId %in% duplicate_ids, ]
  }

  rownames(DT_dat)=DT_dat$SampleId
  condition=DT_dat %>%
    dplyr::select(-matches("seq\\.", ignore.case = TRUE))%>%
    dplyr::rename(SampleID = SampleId)
  condition$SampleID <- as.character(condition$SampleID)
  
  DT_dat = DT_dat[, grep('seq\\.', colnames(DT_dat))] %>%
    t() %>%
    as.data.frame()

  number_samples <- length(names(DT_dat))

  Buffer_mean <- DT %>%
    dplyr::filter(SampleType == "Buffer") %>%
    dplyr::select(matches("seq\\.", ignore.case = TRUE)) %>%
    dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    t() %>%
    as.data.frame()%>%
    {colnames(.) <- 'Buffer'; .}

  Calibrator_mean <- DT %>%
    dplyr::filter(SampleType == "Calibrator") %>%
    dplyr::select(matches("seq\\.", ignore.case = TRUE)) %>%
    dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    t() %>%
    as.data.frame()%>%
    {colnames(.) <- 'Calibrator'; .}

  DT_combined <- cbind(Buffer_mean, Calibrator_mean, DT_dat)
  DT_out=merge(anno[,grep('AptName|UniProt|EntrezGeneSymbol|TargetFullName|Organism|Type',colnames(anno))], DT_combined,by.x='AptName',by.y=0)
  DT_out=DT_out%>%
    dplyr::rename(Protein_Group= UniProt)%>%
    dplyr::rename(Genes= EntrezGeneSymbol)
  return(list(data = DT_out, condition = condition, number_samples = number_samples))
}

Buffer_filter=function(DT){
  DT=as.data.frame(DT)
  DT_filter <- DT %>%
    dplyr::mutate(across(
      .cols = -c(Protein_Group, Genes, Buffer,Calibrator),  # Exclude PG_group, genes, and Buffer
      .fns = ~ ifelse(. < Buffer, NA, .)  # Apply the condition
    ))
  return(DT_filter)
}


#' Filter SummarizedExperiment based on Limit of Detection (LOD)
#'
#' This function filters a SummarizedExperiment object by setting values to NA
#' if they fall below a specified Limit of Detection (LOD). It searches for the
#' LOD values first in the rowData, and then in the assay columns (samples).
#'
#' @param se A SummarizedExperiment object (e.g., proteomics data).
#' @param lod_col A character string specifying the name of the LOD column.
#'                Defaults to "Buffer".
#'
#' @return The modified SummarizedExperiment object with values < LOD set to NA.
#' @export
lod_filter <- function(se, lod_col = "Buffer") {
  
  # 1. Input Validation
  if (!is(se, "SummarizedExperiment")) {
    stop("Input 'se' must be a SummarizedExperiment object.")
  }
  
  lod_values <- NULL
  source_found <- "none"
  
  # 2. Search in rowData (metadata for rows/proteins)
  if (lod_col %in% colnames(rowData(se))) {
    lod_values <- rowData(se)[[lod_col]]
    source_found <- "rowData"
  } 
  # 3. Search in Assay Columns (specific sample acting as LOD, e.g., Buffer)
  else if (lod_col %in% colnames(se)) {
    lod_values <- assay(se, 1)[, lod_col]
    source_found <- "assay"
  } 
  else {
    stop(paste("Could not find", lod_col, "in rowData or assay columns of the SummarizedExperiment."))
  }
  
  # Check if lod_values is numeric
  if (!is.numeric(lod_values)) {
    stop(paste("The values in", lod_col, "are not numeric and cannot be used for filtering."))
  }
  
  # 4. Apply Filter to all assays
  # We iterate through all assays in the object (e.g., counts, log-intensities)
  assay_names <- assayNames(se)
  if (is.null(assay_names)) assay_names <- paste0("assay_", seq_along(assays(se)))
  
  mat <- assay(se)
  mask <- mat < lod_values
  mask[is.na(mask)] <- FALSE
  mat[mask] <- NA
  assay(se) <- mat
  
  return(se)
}

