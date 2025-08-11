library(shiny)
options(shiny.maxRequestSize=30*1024^2)
source("global.R")
source("ui.R")
source("server.R")
source("helpers.R")

shinyApp(ui = ui, server = server)
