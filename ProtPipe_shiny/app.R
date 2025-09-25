library(shiny)
source("global.R")
source("ui.R")
source("server.R")
source("helpers.R")
options(shiny.maxRequestSize=5000 * 1024^2)
set.seed(100)

shinyApp(ui = ui, server = server)
