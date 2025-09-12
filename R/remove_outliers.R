#' Filter Outlier Samples Based on Protein Counts
#'
#' @description
#' This function identifies and removes entire samples that are considered outliers
#' based on the total number of non-missing proteins identified within them.
#'
#' The method calculates the mean and standard deviation of protein counts across
#' all samples. A sample is flagged as an outlier if its total protein count
#' falls outside the range defined by a specified number of standard deviations
#' (`sds`) from the mean.
#'
#' @param object A `ProtData` object.
#' @param sds The number of standard deviations from the mean protein count to use
#'   as the outlier threshold. Defaults to 3.
#'
#' @return A `ProtData` object with the outlier samples removed from the `data`
#'   and `condition` slots.
#'
#' @export
#'
#' @examples
#' # Create data where one sample has a very low number of identified proteins
#' raw_data <- data.frame(
#'   Gene = c("GENEA", "GENEB", "GENEC", "GENED"),
#'   SampleA = c(10, 11, 12, 13),
#'   SampleB = c(12, 13, 14, 15),
#'   SampleC = c(11, 12, 13, 14),
#'   SampleD_outlier = c(NA, NA, NA, 5) # Only one protein found
#' )
#'
#' pd_obj <- create_protdata(dat = raw_data)
#' cat("Dimensions before removing outliers:", dim(pd_obj@data), "\n")
#'
#' # With sds=1, the outlier sample should be identified and removed.
#' filtered_obj <- remove_outliers(pd_obj, sds = 1)
#' cat("Dimensions after removing outliers:", dim(filtered_obj@data), "\n")
#' print(colnames(filtered_obj@data))
#'
setGeneric("remove_outliers", function(object, sds = 3) standardGeneric("remove_outliers"))

setMethod("remove_outliers",
          "ProtData",
          function(object, sds = 3){
            N_values <- colSums(!is.na(getData(object)))
            pgcounts <- data.frame(Sample = names(N_values), N = N_values)

            dat <- object@data
            cond <- object@condition

            stdev <- sd(pgcounts[,'N'])
            mean_count <- mean(pgcounts[,'N'])
            min_protein_groups <- floor(mean_count - (sds * stdev))
            max_protein_groups <- ceiling(mean_count + (sds * stdev))
            cat(paste0('INFO: Tolerating protein group counts in the range [', min_protein_groups,',',max_protein_groups,']\n'))
            low_count_samples <- as.character(pgcounts[pgcounts$N < min_protein_groups, 'Sample'])
            if(length(low_count_samples)>0){
              print("removing the following samples with low protein counts:")
              print(low_count_samples)
            }
            high_count_samples <- as.character(pgcounts[pgcounts$N > max_protein_groups, 'Sample'])
            if(length(high_count_samples)>0){
              print("removing the following samples with high protein counts:")
              print(high_count_samples)
            }

            outliers <- c(low_count_samples, high_count_samples)
            if (length(outliers) >0){
              dat <- dat[,!(colnames(dat) %in% outliers),drop = FALSE]
              cond <- cond[!rownames(cond) %in% outliers, ,drop = FALSE]

              object@data <- dat
              object@condition <- cond
            }
            return(object)
          }
)
