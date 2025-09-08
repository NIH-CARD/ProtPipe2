options(shiny.maxRequestSize=5000 * 1024^2)
library(shiny)
source("global.R")
source("ui.R")
source("server.R")
source("helpers.R")

shinyApp(ui = ui, server = server)
