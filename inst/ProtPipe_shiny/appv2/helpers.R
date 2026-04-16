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
      output$file_ui <- renderUI({
        fileInput(ns("file"), label)
      })
    })

    file
  })
}

add_zip_plot <- function(plot, plotname, subfolder, zip_workspace, zip_file) {
  file_path <- file.path(zip_workspace, subfolder, plotname)
  ggplot2::ggsave(file_path, plot = plot, device = "pdf")
}

add_zip_tabular <- function(data, filename, subfolder, zip_workspace, zip_file) {
  file_path <- file.path(zip_workspace, subfolder, filename)
  data.table::fwrite(data, file_path, sep = "\t")
}

detect_olink_npx <- function(file_path) {
  app_require_packages("OlinkAnalyze", feature = "Olink file detection")

  tryCatch({
    data <- suppressWarnings(OlinkAnalyze::read_NPX(filename = file_path))
    is.data.frame(data) && nrow(data) > 0
  }, error = function(e) {
    FALSE
  })
}
