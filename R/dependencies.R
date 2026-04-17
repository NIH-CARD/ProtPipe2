protpipe_require_packages <- function(packages, feature = "This feature") {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) == 0) {
    return(invisible(TRUE))
  }

  stop(
    sprintf(
      "%s requires the following package%s to be installed: %s",
      feature,
      if (length(missing) == 1) "" else "s",
      paste(missing, collapse = ", ")
    ),
    call. = FALSE
  )
}

protpipe_require_app_packages <- function() {
  protpipe_require_packages(
    c(
      "shiny",
      "shinyjs",
      "bslib",
      "markdown"
    ),
    feature = "The ProtPipe Shiny app"
  )
}
