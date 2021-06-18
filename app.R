library(shiny)
library(shinydashboard)

header <- dashboardHeader()
sidebar <- dashboardSidebar()
body <- dashboardBody()

ui <- dashboardPage(
    header = header,
    sidebar = sidebar,
    body = body
)

server <- function(input, output){
    
}

shinyApp(ui, server)