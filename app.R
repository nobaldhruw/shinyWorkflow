library(shiny)
library(shinydashboard)

# Give a title to your app here
header <- dashboardHeader(
    title = "Differential Analysis"
)

# Put all the menu items in the sidebar
sidebar <- dashboardSidebar(
    sidebarMenu(
        menuItem("Histogram", tabName = "tab_histogram")
    )
)

# Main panel where all the action happens
body <- dashboardBody(
    tabItems(
        tabItem(
            tabName = "tab_histogram",
            fluidRow(
                box(
                    title = 'Inputs', 
                    status = "primary", 
                    width = 4, 
                    solidHeader = TRUE,
                    sliderInput("num", "Choose a number", value = 50, min = 10, max = 100, step = 10)
                ),
                box(
                    title = 'Histogram', 
                    status = "primary", 
                    width = 8, 
                    solidHeader = TRUE,
                    plotOutput("hist")
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
    output$hist <- renderPlot({
        hist(rnorm(input$num))
    })
}

shinyApp(ui, server)