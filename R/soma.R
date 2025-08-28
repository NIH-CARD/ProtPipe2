#' @importFrom magrittr %>%

#' A protData contructor for SomaScan Data
#'
#' @param adat
#'
#' @return
#' @export
#'
#' @examples
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
  if(!is.null(condition)){
    soma_condition <- soma_condition %>%
      left_join(condition, by = "SampleId")
  }
  return(create_protdata(dat, intensity_cols = c(5:length(colnames(dat))), soma_condition, method = "SomaScan"))
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

