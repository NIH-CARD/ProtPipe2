#'
#' This function takes a data frame of proteomics data and its corresponding
#' sample metadata to construct a SummarizedExperiment object. It handles
#' detection of intensity columns, validation, and synchronization of metadata.
#'
#' @param data A data frame containing both protein metadata and intensity values.
#' @param sample_metadata A data frame for sample metadata. It must contain a
#'   'SampleID' column that matches the intensity column headers in 'data'. If
#'   NULL (the default), a basic metadata table is generated from the column names.
#' @param intensity_cols An optional character or numeric vector specifying which
#'   columns in 'data' are the intensity/abundance columns. If NULL, the function

#'   will attempt to autodetect them as all numeric columns.
#' @param creation_method A character string to log how the data was created
#'   (e.g., "MaxQuant_LFQ"). Stored in the object's metadata.
#'
#' @return A \code{SummarizedExperiment} object.
#' @export
#'
#' @examples
#' # --- Sample Data ---
#' test_data <- data.frame(
#'   ProteinID = c("P02768", "P01023", "P10636"),
#'   Gene = c("ALB", "A2M", "CALM1"),
#'   Sample_A_1 = c(10.1, 12.5, 13.2),
#'   Sample_A_2 = c(10.5, 12.8, 13.9),
#'   Sample_B_1 = c(15.2, 18.1, 19.0),
#'   Sample_B_2 = c(15.8, 18.5, 19.9)
#' )
#'
#' sample_info <- data.frame(
#'   SampleID = c("Sample_A_1", "Sample_A_2", "Sample_B_1", "Sample_B_2"),
#'   Condition = c("A", "A", "B", "B")
#' )
#'
#' # --- Function Call ---
#' se <- create_se(
#'   data = test_data,
#'   sample_metadata = sample_info
#' )
#'
#' # --- Inspect the new object ---
#' print(se)
#' assayNames(se)
#' head(rowData(se))
#' colData(se)
#' metadata(se)
#'
create_se <- function(data, sample_metadata = NULL, intensity_cols = NULL, creation_method = "Unknown") {

  # --- 1. Initial Checks and Setup ---
  if (!is.data.frame(data)) {
    stop("The 'data' argument must be a data frame.")
  }

  data <- convert_numeric_cols(data)
  colnames(data) <- trim_names(colnames(data))

  # Autodetect intensity columns if not provided
  if (is.null(intensity_cols)) {
    message("`intensity_cols` not provided. Detecting numeric columns as intensity data.")
    intensity_cols <- which(sapply(data, is.numeric))
    if (length(intensity_cols) == 0) {
      stop("No numeric intensity columns could be detected.")
    }
  }

  # --- 2. Separate Data into Assay, RowData, and ColData Components ---

  # The assay matrix (must be numeric)
  assay_data <- as.matrix(data[, intensity_cols])

  # The row metadata (protein metadata)
  # Ensure it's a DataFrame for Bioconductor consistency
  row_data <- S4Vectors::DataFrame(data[, -intensity_cols, drop = FALSE])

  # --- 3. Process and Synchronize Sample Metadata (colData) ---
  if (!is.null(sample_metadata)) {
    col_data <- as.data.frame(sample_metadata)

    if (!"SampleID" %in% colnames(col_data)) {
      stop("Error: The sample_metadata file is missing the required 'SampleID' column.")
    }
    if (any(duplicated(col_data$SampleID))) {
      stop("Error: The 'SampleID' column contains duplicate values.")
    }

    # Set rownames from the SampleID column
    rownames(col_data) <- trim_names(as.character(col_data$SampleID))
    col_data$SampleID <- NULL

    # Synchronize rownames of col_data with colnames of assay_data
    assay_colnames <- colnames(assay_data)

    # Check for mismatches
    if (!setequal(rownames(col_data), assay_colnames)) {
      warning("SampleIDs in metadata do not perfectly match intensity column names.")

      # Drop rows of col_data that are not in data
      matching_rows <- intersect(rownames(col_data), colnames(assay_data))
      col_data <- col_data[matching_rows, ]

      dropped_rows <- setdiff(rownames(col_data), matching_rows)
      warning("Dropped rows:\n")
      warning(dropped_rows)
    }
    # add additional rows to col_data if they exist in data
    if (ncol(assay_data) > nrow(col_data)) {
      # Find missing columns
      missing_cols <- setdiff(colnames(assay_data), rownames(col_data))

      # Create rows with NA for the missing columns and add them to 'col_data'
      missing_rows <- data.frame(matrix(NA, nrow = length(missing_cols), ncol = ncol(col_data)))
      colnames(missing_rows) <- names(col_data)
      rownames(missing_rows) <- missing_cols

      # Add missing rows to 'col_data'
      col_data <- rbind(col_data, missing_rows)
    }

      #assay_data <- assay_data[, shared_samples, drop = FALSE]
      #col_data <- col_data[shared_samples, , drop = FALSE]

    # Reorder col_data to perfectly match the order of assay_data
    col_data <- col_data[colnames(assay_data), , drop = FALSE]

  } else {
    # If no sample_metadata is provided, create a default one
    warning("`sample_metadata` not provided. Generating a basic version from column names.")
    col_data <- data.frame(
      row.names = colnames(assay_data),
      base_condition = gsub("_\\d+$", "", colnames(assay_data))
    )
  }

  # ensure proper sample names are used
  colnames(assay_data) <- make.names(colnames(assay_data))
  colnames(col_data) <- make.names(colnames(col_data))
  rownames(col_data) <- make.names(rownames(col_data))

  # --- 4. Construct the SummarizedExperiment Object ---

  # IMPORTANT: The commented-out dplyr block from your original function for
  # filtering proteins should NOT be in a constructor. That is a downstream
  # *processing step*. This function should only build the object.

  se <- SummarizedExperiment(
    # Assays must be in a named list. We'll call the primary data "intensities".
    assays = list(intensities = assay_data),

    # Protein metadata
    rowData = row_data,

    # Sample metadata
    colData = col_data,

    # Experiment-level metadata, a perfect place for your processing log
    metadata = list(
      creation_method = creation_method,
      processing_log = list() # Initialize the log
    )
  )

  return(se)
}

#' Identify Numeric Columns by Index
#'
#' @description
#' A helper function that scans a data frame and returns the integer indices
#' of the columns that contain numeric data.
#'
#' This is primarily used internally by data loading functions to automatically
#' distinguish quantitative intensity columns from metadata columns.
#'
#' @param df A `data.frame` to be scanned.
#'
#' @return An integer vector containing the column indices of all numeric
#'   columns found in the input `data.frame`.
#'
#' @keywords internal
#'
#' @export
detect_intensity_cols <- function(df) {
  which(sapply(df, is.numeric))
}

#' Check if a Processing Step has been Applied
#'
#' Scans an object's processing history to determine if a step with a specific
#' name has already been performed. This is useful for preventing redundant
#' operations or for enforcing a required order of preprocessing steps.
#'
#' @param object An S4 object containing a `processing` slot (list).
#' @param step_name A character string specifying the unique name of the step to
#'   check for (e.g., `"log2_transform"`).
#'
#' @return A logical value: `TRUE` if a log entry with the given `step_name`
#'   exists in the object's history, `FALSE` otherwise.
#'
#' @export
has_step <- function(object, step_name) {
  if (length(metadata(object)$processing_log) == 0) return(FALSE)

  # Get the names of all steps performed so far
  executed_steps <- sapply(metadata(object)$processing_log, function(x) x$name)

  return(step_name %in% executed_steps)
}

######## HELPER FUNCTIONS ####################################################################

#total proteomics############
standardize_format <- function(DT.original) {
  # Accepts an input protein group intensity data.table, whether spectronaut or DIA-NN format,
  # and restructures into one consistent style for downstream processing
  DT <- DT.original
  if("Protein.Ids" %in% colnames(DT) | "Protein.Group" %in% colnames(DT)) {
    # print("DIAnn input")
    # DT[, 'Protein.Ids' := NULL]
    # DT[, 'Protein.Names' := NULL]
    # DT[, 'First.Protein.Description' := NULL]
    data.table::setnames(DT, 'Protein.Group', 'Protein_Group')
  }
  else if('EG.PrecursorId' %in% colnames(DT)) {
    print("Spectronaut input")
    data.table::setnames(DT, 'EG.PrecursorId', 'Peptide_Sequence')
    data.table::setnames(DT, 'PG.Genes', 'Genes')
    DT=as.data.frame(DT)
    # Use only Protein_Group and Genes
    dplyr::select=c('Peptide_Sequence','Genes',grep('raw',colnames(DT),value = T))
    DT=DT[, dplyr::select]
    #as number
    DT[,grep('raw',colnames(DT))]=as.data.frame(apply(DT[,grep('raw',colnames(DT))],2,as.numeric))
    DT=data.table(DT)
  }
  else if('PG.ProteinGroups' %in% colnames(DT)) {
    print("Spectronaut input")
    data.table::setnames(DT, 'PG.ProteinGroups', 'Protein_Group')
    data.table::setnames(DT, 'PG.Genes', 'Genes')
    DT=as.data.frame(DT)
    # Use only Protein_Group and Genes
    dplyr::select=c('Protein_Group','Genes',grep('PG.Quantity',colnames(DT),value = T))
    DT=DT[, dplyr::select]
    #as number
    DT[,grep('PG.Quantity',colnames(DT))]=as.data.frame(apply(DT[,grep('PG.Quantity',colnames(DT))],2,as.numeric))
    DT=data.table(DT)
  }
  else if('Peptide Sequence' %in% colnames(DT)) {
    print("FragPipe input")
    data.table::setnames(DT, 'Gene', 'Genes')
    colnames(DT)=gsub("\\s", "_",colnames(DT))
    # Use only Protein_Group and Genes
    DT=as.data.frame(DT)
    dplyr::select=c('Peptide_Sequence','Genes',grep('[0-9]_Intensity',colnames(DT),value = T))
    DT=DT[, dplyr::select]
    DT=data.table(DT)
  }

  # Remove leading directories for sample names
  # e.g. /path/to/sample1.mzML -> sample1.mzML
  data.table::setnames(DT, basename(colnames(DT)))

  # Remove trailing file extensions
  extensions <- '.mzML$|.mzml$|.RAW$|.raw$|.dia$|.DIA$|_Intensity'
  extension_samplenames <-  colnames(DT)[data.table::`%like%`(colnames(DT), extensions)]
  trimmed_samplenames <- gsub(extensions, '', extension_samplenames)
  data.table::setnames(DT, extension_samplenames, trimmed_samplenames)
  return(DT[])
}

trim_names <- function(names) {
  colnames_out <- gsub(pattern="\\[.*\\] ", replacement='', x=names)   # trim leading [N]
  colnames_out <- gsub(pattern="\\..*\\.PG\\.Quantity|\\.PG\\.Quantity|\\..*Quantity.*", replacement='', x=colnames_out)   # remove suffix
  # Remove everything before the last "/" and remove extensions like .raw or .mzml
  colnames_out <- gsub(pattern=".*/", replacement='', x=colnames_out)
  colnames_out <- gsub(pattern="\\.(raw|mzML)$", replacement='', x=colnames_out)
  return(colnames_out)
}

#' Safely Convert Character Columns to Numeric Type
#'
#' @description
#' This function inspects each column of a data frame. If a character column
#' contains only values that can be interpreted as numbers (including integers,
#' decimals, scientific notation, and the strings "NA" or "NaN"), it converts
#' that entire column to the numeric type. Columns containing non-numeric
#' text are left unchanged.
#'
#' @param df A data frame whose character columns will be checked for potential
#'   conversion to numeric.
#'
#' @return The input data frame, with any eligible character columns converted
#'   to the numeric class.
#'
#' @export
#'
convert_numeric_cols <- function(df) {

  # Regex for a valid number string (handles integers, decimals, scientific notation)
  numeric_regex <- "^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$"

  df[] <- lapply(df, function(col) {
    # Only consider non-numeric (character) columns for conversion
    if (is.character(col)) {

      # Remove true R NA values and empty strings for the check.
      # We want to evaluate the actual character strings present.
      values_to_check <- col[!is.na(col) & col != ""]

      # If there are no actual strings to check, no need to convert
      if (length(values_to_check) == 0) {
        return(col)
      }

      # THE CRUCIAL CHECK:
      # Every string must either match the numeric regex OR be "NA" or "NaN".
      # We use toupper() to make the NA/NaN check case-insensitive.
      is_valid_entry <- grepl(numeric_regex, values_to_check) | toupper(values_to_check) %in% c("NA", "NAN")

      if (all(is_valid_entry)) {
        # If all entries are valid, it's safe to convert the entire column.
        # as.numeric() correctly handles the strings "NA" and "NaN".
        return(as.numeric(col))
      }
    }

    # If the column is not character or fails the check, return it unmodified.
    return(col)
  })

  return(df)
}

#create long data table
melt_intensity_table <- function(DT) {
  # Converts intensity data.table to long format
  # info_cols <- c('Protein_Group', 'Genes', 'First_Protein_Description')
  DT.long <- reshape2::melt(DT,
                  measure.vars=names(DT)[sapply(DT, function(x) all(is.numeric(x)))],
                  variable.name='Sample',
                  value.name='Intensity')
  DT.long=data.table::data.table(DT.long)
  DT.long=DT.long[DT.long$Intensity>0,]
  return(DT.long)
}
