library(shiny)
fileUploadUI <- function(id, label = "Upload File") {
  ns <- NS(id)
  tagList(
    fileInput(ns("file"), label),
    actionButton(ns("clear"), "Remove file")
  )
}

fileUploadServer <- function(id, label = "Upload File") {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    file <- reactiveVal(NULL)

    output$file_ui <- renderUI({
      fileInput(ns("file"), label)
    })

    observeEvent(input$file, {
      file(input$file)
    })

    observeEvent(input$clear, {
      file(NULL)
      # Force UI to re-render, fully resetting the file input
      output$file_ui <- renderUI({
        fileInput(ns("file"), label)
      })
    })

    return(file)
  })
}

add_zip_plot <- function(plot, plotname, subfolder, zip_workspace, zip_file){
  file_path <- file.path(zip_workspace, subfolder, plotname)
  ggsave(file_path, plot=plot, device = "pdf")
  # rel_path <- file.path(subfolder, plotname)
  # zip::zip_append(
  #   zipfile = file.path(zip_workspace, zip_file),
  #   files = rel_path,
  #   root = zip_workspace
  # )
  # #file.remove(file_path)
}

add_zip_tabular <- function(data, filename, subfolder, zip_workspace, zip_file){
  file_path <- file.path(zip_workspace, subfolder, filename)
  data.table::fwrite(data, file_path, sep = "\t")
  # rel_path <- file.path(subfolder, filename)
  # zip::zip_append(
  #   zipfile = file.path(zip_workspace, zip_file),
  #   files = rel_path,
  #   root = zip_workspace
  # )
  # file.remove(file_path)
}

detect_olink_npx <- function(file_path) {

  # 2. Validate that the file exists
  if (!file.exists(file_path)) {
    warning("File does not exist: ", file_path)
    return(FALSE)
  }

  # 3. Read only the header based on the file extension
  ext <- tolower(tools::file_ext(file_path))

  header_df <- tryCatch({
    if (ext %in% c("xls", "xlsx")) {
      # For Excel files, use readxl to read zero rows
      readxl::read_excel(file_path, n_max = 0)
    } else if (ext %in% c("csv", "tsv", "txt")) {
      # For text files, use data.table::fread to read zero rows
      data.table::fread(file_path, nrows = 0)
    } else {
      warning("Unsupported file extension: ", ext)
      return(NULL)
    }
  }, error = function(e) {
    warning("Failed to read file header. Error: ", e$message)
    return(NULL)
  })

  # 4. If reading the header failed, it's not a valid file
  if (is.null(header_df)) {
    return(FALSE)
  }

  # 5. Check for the Olink signature columns
  file_colnames <- colnames(header_df)
  olink_signature_cols <- c("OlinkID", "Panel", "NPX")

  return(all(olink_signature_cols %in% file_colnames))
}
