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
  ggplot2::ggsave(file_path, plot=plot, device = "pdf")
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

library(OlinkAnalyze)

detect_olink_npx <- function(file_path) {
  # Attempt to read the file using the official Olink loader
  # We suppress warnings to keep the output clean if the read fails
  tryCatch({
    # Try reading just a small subset if possible, but OlinkAnalyze usually needs to
    # parse the whole format to validate. The key is that read_NPX handles zip/csv/txt automatically.
    data <- suppressWarnings(OlinkAnalyze::read_NPX(filename = file_path))
    
    # Check if the result is a valid tibble/data.frame and has content
    return(is.data.frame(data) && nrow(data) > 0)
    
  }, error = function(e) {
    # If read_NPX throws an error, it is not a valid Olink file
    return(FALSE)
  })
}
