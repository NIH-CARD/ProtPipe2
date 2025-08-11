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
