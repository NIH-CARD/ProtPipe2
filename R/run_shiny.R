#' Run the protpipe Shiny App
#'
#' This function launches the protpipe_shiny application stored within the
#' package. It handles locating the app directory and ensuring the
#' necessary packages are loaded.
#'
#' @export
#' @return A Shiny application object (run in a separate R session/browser).
run_protpipe_shiny <- function() {
  # Find the path to the 'protpipe_shiny' directory within the package
  appDir <- system.file("ProtPipe_shiny", package = "ProtPipe")

  if (appDir == "") {
    stop("Could not find the 'protpipe_shiny' directory. Try re-installing the protpipe package.")
  }

  # Ensure the shiny package is loaded (it should be a dependency)
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The 'shiny' package is required to run this application. Please install it using install.packages('shiny').")
  }

  # Launch the application
  shiny::runApp(appDir, display.mode = "normal")
}
