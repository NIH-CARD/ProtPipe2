#' @importFrom magrittr %>%

#' Remove Specific Samples from a ProtData Object
#'
#' @description
#' This function removes one or more specified samples from all relevant slots
#' of a `ProtData` object, including the quantitative data (`data` slot) and
#' the sample condition metadata (`condition` slot).
#'
#' @param object A `ProtData` object.
#' @param samples A character vector containing the exact names of the samples to
#'   be removed.
#'
#' @return A `ProtData` object with the specified samples removed.
#'
#' @export
#'
#' @examples
#' # Create a sample ProtData object with four samples
#' raw_data <- data.frame(
#'   Gene = c("GENEA", "GENEB", "GENEC"),
#'   SampleA = c(10, 11, 12),
#'   SampleB = c(12, 13, 14),
#'   SampleC = c(11, 12, 13),
#'   SampleD = c(14, 15, 16)
#' )
#' pd_obj <- create_protdata(dat = raw_data)
#' cat("Sample names before removal:\n")
#' print(colnames(pd_obj@data))
#'
#' # Remove two of the samples by name
#' filtered_obj <- remove_sample(pd_obj, samples = c("SampleA", "SampleC"))
#'
#' cat("\nSample names after removal:\n")
#' print(colnames(filtered_obj@data))
#'
setGeneric("remove_sample", function(object, samples) standardGeneric("remove_sample"))

setMethod("remove_sample",
          "ProtData",
          function(object, samples){
            dat <- object@data
            cond <- object@condition

            dat <- dat[,!(colnames(dat) %in% samples),drop = FALSE]
            cond <- cond[!rownames(cond) %in% samples, ,drop = FALSE]

            object@data <- dat
            object@condition <- cond
            return(object)
          }
)
