options(shiny.sanitize.errors = TRUE)
library(shiny)
library(shinyjs)
library(bslib)

app_require_packages <- function(packages, feature = "This app feature") {
  ProtPipe:::protpipe_require_packages(packages, feature = feature)
}

organism_map <- list(
  "Human" = list(orgdb_package = "org.Hs.eg.db", kegg = "hsa"),
  "Mouse" = list(orgdb_package = "org.Mm.eg.db", kegg = "mmu"),
  "Rat" = list(orgdb_package = "org.Rn.eg.db", kegg = "rno"),
  "Fly" = list(orgdb_package = "org.Dm.eg.db", kegg = "dme"),
  "Nematode" = list(orgdb_package = "org.Ce.eg.db", kegg = "cel")
)
