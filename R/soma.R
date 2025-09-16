#' Create a ProtData Object from SomaScan Data
#'
#' @description
#' Processes SomaScan data from an ADAT object, optionally performs filtering
#' (e.g., on buffer controls), and formats the result into a ProtData object
#' for downstream analysis.
#'
#' @param adat A SomaScan ADAT object, typically loaded using the SomaDataIO package.
#' @param condition An optional data frame containing sample metadata. This is passed to
#'   the downstream create_protdata function for integration.
#' @param filter A logical value (TRUE/FALSE) indicating whether to perform
#'   filtering on the data, including for buffer controls. Defaults to TRUE.
#'
#' @return A ProtData object formatted for use with other package functions.
#'
#' @export
#'
#' @importFrom dplyr rename
#' @importFrom magrittr %>%
#'
create_protdata_from_soma <- function(adat, condition = NULL, filter = TRUE) {

  if(filter){
    soma_out <- soma_all_output(adat)
    dat <- soma_out$data
    soma_condition <- soma_out$condition
    dat <- Buffer_filter(dat)
  }else{
    soma_out <- soma_sample_out(adat)
    dat <- soma_out$data
    soma_condition <- soma_out$condition
  }
  soma_condition <- soma_condition %>%
    dplyr::rename(SampleID = SampleId)
  # if(!is.null(condition)){
  #   soma_condition <- soma_condition %>%
  #     dplyr::left_join(condition, by = "SampleID")
  # }
  return(create_protdata(dat, condition = condition, method = "SomaScan"))
}


soma_sample_out=function(DT){
  anno <- SomaDataIO::getAnalyteInfo(DT)%>%
    dplyr::filter(Organism == "Human") %>%
    dplyr::filter(Type == "Protein")
  DT=as.data.frame(DT)
  DT_dat=DT%>%
    dplyr::filter(grepl("Sample", SampleType, ignore.case = TRUE))

  #check for duplicated SampleID
  duplicate_ids <- DT_dat$SampleId[duplicated(DT_dat$SampleId) | duplicated(DT_dat$SampleId, fromLast = TRUE)]
  if(length(duplicate_ids>0)){
    cat(paste0("removing duplicates: ", paste(duplicate_ids, collapse = ", ")))
    DT_dat <- DT_dat[!DT_dat$SampleId %in% duplicate_ids, ]
  }


  rownames(DT_dat)=DT_dat$SampleId
  condition=DT_dat %>%
    dplyr::select(-matches("seq\\.", ignore.case = TRUE))
  DT_dat=DT_dat%>%
    dplyr::select(matches("seq\\.", ignore.case = TRUE))%>%
    t()
  DT_dat=merge(anno[,grep('AptName|UniProt|EntrezGeneSymbol|TargetFullName',colnames(anno))], DT_dat,by.x='AptName',by.y=0,all.x=T)
  DT_dat=DT_dat %>%
    dplyr::filter(UniProt != "")%>%
    dplyr::filter(EntrezGeneSymbol != "") %>%
    dplyr::rename(Protein_Group= UniProt)%>%
    dplyr::rename(Genes= EntrezGeneSymbol)
  return(list(data = DT_dat, condition = condition))
}

##format data
soma_all_output=function(DT){
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
    dplyr::select(-matches("seq\\.", ignore.case = TRUE))
  DT_dat = DT_dat[, grep('seq\\.', colnames(DT_dat))] %>%
    t() %>%
    as.data.frame()

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
  return(list(data = DT_out, condition = condition))
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

