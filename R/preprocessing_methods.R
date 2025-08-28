#' Define the get_overlap generic method
#'
#' @param object The object to operate on.
#' @param condition_name The name of the condition to find overlap across.
#' @export
setGeneric("get_overlap",
           def = function(object, condition_name) {
             standardGeneric("get_overlap")
           }
)

#' get_overlap method for protdata class
#'
#' Filters the protdata object to retain only proteins that are present (non-NA)
#' in every sample within each unique group of the specified condition.
#'
#' @param object A protdata object.
#' @param condition_name A character string specifying the column name in the
#' condition' slot to group samples by.
#'
#' @return A new, modified protdata object containing only the overlapping proteins.
setMethod("get_overlap",
          signature(object = "ProtData", condition_name = "character"),
          function(object, condition_name) {

            # --- Input Validation ---
            if (length(condition_name) != 1) {
              stop("'condition_name' must be a single character string.")
            }
            if (!(condition_name %in% names(object@condition))) {
              stop(paste0("'", condition_name, "' is not a valid column in the 'condition' slot."))
            }

            cat("Starting overlap analysis for condition: '", condition_name, "'\n", sep = "")

            # --- Extract data from slots ---
            cond_df <- object@condition
            data_df <- object@data

            # Get the unique groups from the specified condition column
            groups <- unique(cond_df[[condition_name]])

            # A list to store the names of present proteins for each group
            present_proteins_per_group <- list()

            # --- Logic to find present proteins for each group ---
            for (group in groups) {
              # Get the sample names belonging to the current group
              samples_in_group <- rownames(cond_df)[cond_df[[condition_name]] == group]

              # Subset the data to include only these samples
              # Use drop=FALSE to ensure it remains a data.frame even with one sample
              data_subset <- data_df[, samples_in_group, drop = FALSE]

              # A protein is "present" if it has no NA values across all samples in this group
              # This is a strict definition of presence.
              is_present <- rowSums(!is.na(data_subset)) > 0

              # Get the names of the proteins that are present
              present_protein_names <- rownames(data_subset)[is_present]

              present_proteins_per_group[[as.character(group)]] <- present_protein_names

              cat("  - Group '", group, "': Found ", length(present_protein_names), " present proteins.\n", sep = "")
            }

            # --- Find the intersection of all sets of proteins ---
            if (length(present_proteins_per_group) == 0) {
              stop("No groups found for the specified condition.")
            }

            # Use Reduce with intersect to find common proteins across all lists
            overlapping_proteins <- Reduce(intersect, present_proteins_per_group)

            cat("Found", length(overlapping_proteins), "proteins present across all groups.\n")

            if (length(overlapping_proteins) == 0) {
              warning("No overlapping proteins found. Returning an empty object.")
            }

            # --- Modify the object slots based on the overlap ---
            # Note: In R, we return a modified copy rather than modifying in-place.
            # The user will typically assign the result back to the original variable.

            object@prot_meta <- object@prot_meta[overlapping_proteins, , drop = FALSE]
            object@data <- object@data[overlapping_proteins, , drop = FALSE]

            cat("Object has been updated.\n")

            return(object)
          }
)
