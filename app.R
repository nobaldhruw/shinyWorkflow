library(shiny)
library(shinydashboard)
library(DT)
library(cmapR)

# Define some global variables
gct <- parse.gctx('data/GSE104310.gct')
exp <- gct@mat
colData <- gct@cdesc

# Give a title to your app here
header <- dashboardHeader(
    title = "Differential Analysis"
)

# Put all the menu items in the sidebar
sidebar <- dashboardSidebar(
    sidebarMenu(
        menuItem("Data", tabName = "tab_data")
    )
)

# Main panel where all the action happens
body <- dashboardBody(
    tabItems(
        tabItem(
            tabName = "tab_data",
            fluidRow(
                box(
                    title = 'Expression Data', 
                    status = "primary", 
                    width = 12, 
                    solidHeader = TRUE,
                    dataTableOutput("exp_data_output")
                )
            ),
            fluidRow(
                box(
                    title = 'Sample metadata', 
                    status = "primary", 
                    width = 12, 
                    solidHeader = TRUE,
                    dataTableOutput("colData_output")
                )
            )
        )
    )
)

ui <- dashboardPage(
    header = header,
    sidebar = sidebar,
    body = body
)

server <- function(input, output){
    
    output$exp_data_output <- renderDT(
        exp,
        options = list(pageLength = 5, scrollX=TRUE)
    )
    output$colData_output <- renderDT(
        colData,
        options = list(pageLength = 5, scrollX=TRUE)
    )
}

shinyApp(ui, server)