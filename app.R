library(shiny)
library(shinydashboard)

# Give a title to your app here
header <- dashboardHeader(
    title = "Differential Analysis"
)

# Put all the menu items in the sidebar
sidebar <- dashboardSidebar(
    sidebarMenu(
        menuItem("Data", tabName = "input_data"),
        menuItem("Viz", tabName = "visualization")
    )
)

# Main panel where all the action happens
body <- dashboardBody(
    tabItems(
        tabItem(
            tabName = "input_data",
            h1("Input data"),
            p("This tab will display the data")
        ),
        tabItem(
            tabName = "visualization",
            h1("Visualization"),
            p("This tab will have some graphs/plots")
        )
    )
)

ui <- dashboardPage(
    header = header,
    sidebar = sidebar,
    body = body
)

server <- function(input, output){
    
}

shinyApp(ui, server)