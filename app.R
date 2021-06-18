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
        menuItem("Data", tabName = "tab_data"),
        menuItem("PCA", tabName = "tab_pca")
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
        ),
        tabItem(
            tabName = "tab_pca",
            fluidRow(
                box(
                    title = 'Inputs', 
                    status = "primary", 
                    width = 3, 
                    solidHeader = TRUE,
                    uiOutput("pca_input")
                ),
                box(
                    title = 'Scores plot', 
                    status = "primary", 
                    width = 9, 
                    solidHeader = TRUE,
                    plotOutput("pca_scores_plot")
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
    output$pca_input <- renderUI({
        tagList(
            selectInput("pc_x", "Choose PC on x-axis", choices = 1:5, selected = 1),
            selectInput("pc_y", "Choose PC on y-axis", choices = 1:5, selected = 2),
            selectInput("pca_color", "Choose column to color", choices = unique(colnames(colData)))
        )
    })
    output$pca_scores_plot <- renderPlot({
        PCAObj <- prcomp(as.data.frame(t(exp)))
        PCAScores <- data.frame(PCAObj$x, cohort = as.factor(as.character(colData[,input$pca_color, drop = T])))
        p <- ggplot(PCAScores, aes_string(x = paste0('PC', input$pc_x), 
                                          y = paste0('PC', input$pc_y),
                                          fill = 'cohort'))+
                geom_point(shape=21, size=4, alpha=0.7)
        p
    })
}

shinyApp(ui, server)