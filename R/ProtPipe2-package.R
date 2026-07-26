#' ProtPipe2: proteomics analysis tool
#'
#' Reproducible workflows for downstream proteomics analysis. All exported
#' functions accept and return \code{SummarizedExperiment} objects.
#'
#' @section Imported accessors:
#' Package code calls the \code{SummarizedExperiment} accessors unqualified, so
#' they are imported here rather than repeated across the individual method
#' files. The same applies to the \code{stats}, \code{utils} and \code{dplyr}
#' helpers used throughout the package.
#'
#' Note that \code{data.table} is deliberately *not* imported: doing so would
#' make the package "data.table aware" and change \code{[} semantics package
#' wide, so its functions are called with \code{data.table::} instead.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom SummarizedExperiment assay assay<- assays assayNames rowData colData
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom dplyr across matches
#' @importFrom rlang .data
#' @importFrom methods is
#' @importFrom utils head
#' @importFrom stats setNames as.formula model.matrix median quantile sd var cor cor.test na.omit dist hclust p.adjust
## usethis namespace: end
NULL

# Column names used for non-standard evaluation inside dplyr/data.table verbs.
# Declaring them keeps R CMD check from reporting them as undefined globals.
utils::globalVariables(c(
  ".", "Assay", "Buffer", "CV", "Calibrator", "Condition", "EntrezGeneSymbol",
  "Genes", "Group", "HeatmapCondition", "Intensity", "LOD", "MeanIntensity",
  "NPX", "OlinkID", "Protein", "ProteinID", "Protein_Group", "Protein_Groups",
  "Sample", "SampleA", "SampleB", "SampleID", "SampleId", "SampleType",
  "Spearman", "UMAP1", "UMAP2", "UniProt", "ZScore", "groupcontrol",
  "grouptreatment", "labeltext", "logFC", "median_intensity", "missing_value",
  "org.Hs.eg.db", "p.signif", "pval", "rho", "rn"
))
