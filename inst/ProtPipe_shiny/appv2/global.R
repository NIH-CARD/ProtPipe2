options(shiny.sanitize.errors = TRUE)
options(shiny.maxRequestSize = 5000 * 1024^2)

source("helpers.R")

library(shiny)
app_require_packages <- function(packages, feature = "This app feature") {
  ProtPipe2:::protpipe_require_packages(packages, feature = feature)
}

app_require_packages(
  c("shinyjs", "bslib", "markdown"),
  feature = "The ProtPipe Shiny app"
)

organism_map <- list(
  "Human" = list(orgdb_package = "org.Hs.eg.db", kegg = "hsa"),
  "Mouse" = list(orgdb_package = "org.Mm.eg.db", kegg = "mmu"),
  "Rat" = list(orgdb_package = "org.Rn.eg.db", kegg = "rno"),
  "Fly" = list(orgdb_package = "org.Dm.eg.db", kegg = "dme"),
  "Nematode" = list(orgdb_package = "org.Ce.eg.db", kegg = "cel")
)
