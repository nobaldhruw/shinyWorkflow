library(shiny)
library(shinydashboard)
library(DT)
library(cmapR)

options(shiny.maxRequestSize=1024*1024^2)

# Import helper scripts
source('scripts/helper.R')

# Give a title to your app here
header <- dashboardHeader(
  title = "Differential Analysis"
)

# Put all the menu items in the sidebar
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Data", tabName = "tab_data"),
    menuItem("PCA", tabName = "tab_pca"),
    menuItem("Differential Expression", tabName = "tab_de"),
    menuItem("EnrichR ", tabName = "tab_enrichr"),
    menuItem("X2K", tabName = "tab_x2k")
  )
)

# Main panel where all the action happens
body <- dashboardBody(
  tabItems(
    tabItem(
      tabName = "tab_data",
      fluidRow(
        box(
          title = 'Import data',
          fileInput("inputFile", "Upload GCT file", accept = ".gct")
        )
      ),
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
    ),
    tabItem(
      tabName = "tab_de",
      fluidRow(
        box(
          title = 'Input parameters',
          status = 'primary',
          width = 3,
          solidHeader = TRUE,
          selectInput("de_cohortCol", "Choose cohort column", choices = NULL),
          selectInput("de_cohortA", "Cohort A", choices = NULL, multiple = TRUE),
          selectInput("de_cohortB", "Cohort B", choices = NULL, multiple = TRUE),
          actionButton("de_button", "Go!")
        ),
        box(
          title = 'Differential Expression Result',
          status = 'primary',
          width = 9,
          solidHeader = TRUE,
          dataTableOutput("de_result_table")
        )
      ),
      fluidRow(
        box(
          title = 'Customize plot',
          status = 'primary',
          width = 3,
          solidHeader = TRUE,
          sliderInput("de_abs_logFC_cutoff", "Absolute logFC cutoff", min = 0, max = 10, value = 1),
          numericInput("de_fdr_cutoff", "FDR cutoff", min = 0, max = 1, value = 0.05),
          numericInput("de_ngenes_to_label", "No. of genes to label", min = 1, max = 20, value = 5)
        ),
        box(
          title = 'Volcano plot',
          status = 'primary',
          width = 9,
          solidHeader = TRUE,
          plotOutput("de_volcano_plot")
        )
      )
    ),
    tabItem(
      tabName = "tab_enrichr",
      fluidRow(
        box(
          title = 'Input paramters',
          status = 'primary',
          width = 3,
          solidHeader = TRUE,
          radioButtons("enrichr_gene_direction", "Gene direction", choices=c("up", "down", "both"), selected = "both"),
          sliderInput("enrichr_abs_logFC_cutoff", "Absolute logFC cutoff", min = 0, max = 10, value = 1),
          numericInput("enrichr_fdr_cutoff", "FDR cutoff", min = 0, max = 1, value = 0.05),
          selectInput("enrichr_db", "Choose DB", choices=NULL),
          numericInput("enrichr_ngs", "No. of gene sets", min=1, max=20, value=10)
        ),
        box(
          title = 'Enriched Gene sets',
          status = 'primary',
          width = 9,
          solidHeader = TRUE,
          plotOutput("enrichr_dotplot")
        )
      )
    ),
    tabItem(
      tabName = "tab_x2k",
      fluidRow(
        box(
          title = 'Input paramters',
          status = 'primary',
          width = 3,
          solidHeader = TRUE,
          radioButtons("x2k_gene_direction", "Gene direction", choices=c("up", "down", "both"), selected = "both"),
          sliderInput("x2k_abs_logFC_cutoff", "Absolute logFC cutoff", min = 0, max = 10, value = 1),
          numericInput("x2k_fdr_cutoff", "FDR cutoff", min = 0, max = 1, value = 0.05),
          numericInput("x2k_nbars", "No. of bars to show", min=1, max=20, value=10)
        ),
        box(
          title = 'Enriched TFs/Kinases',
          status = 'primary',
          width = 9,
          solidHeader = TRUE,
          plotOutput("x2k_barplot")
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

server <- function(input, output, session){
  
  rv <- reactiveValues(
    exp=NULL, 
    colData=NULL, 
    pca_scores=NULL, 
    exp.log_transformed=NULL, 
    de_result=NULL,
    enrichr_dbs=get_enrichr_dblist(),
    egs=NULL,
    x2k_result=NULL
  )
  
  observeEvent(input$inputFile, {
    gctObj <- parse_gctx(input$inputFile$datapath)
    rv$exp <- gctObj@mat
    rv$colData <- gctObj@cdesc
    rv$exp.log_transformed <- log2(rv$exp+1)
    
    cohortCols <- colnames(rv$colData)
    updateSelectInput(session, inputId = "de_cohortCol", choices = cohortCols, selected = cohortCols[1])
    updateSelectInput(session, inputId = "enrichr_db", choices = rv$enrichr_dbs, selected = rv$enrichr_dbs[1])
    
  })

  observeEvent(input$de_cohortCol, {
    updateSelectInput(session, inputId = "de_cohortA", choices = unique(rv$colData[, input$de_cohortCol]))
    updateSelectInput(session, inputId = "de_cohortB", choices = unique(rv$colData[, input$de_cohortCol]))
  })
  
  observeEvent(input$pca_button, {
    rv$pca_scores <- compute_pca(rv$exp.log_transformed, rv$colData)
  })
  
  observeEvent(input$de_button, {
    rv$de_result <- diff_exp_limma(rv$exp.log_transformed, rv$colData, input$de_cohortCol, input$de_cohortA, input$de_cohortB, 'BH')
  })
  
  observeEvent(rv$de_result, {
      rv$egs <- get_enriched_gene_sets(rv$de_result, 
                                       gene_direction  = input$enrichr_gene_direction,
                                       log2fc_cutoff   = input$enrichr_abs_logFC_cutoff, 
                                       fdr_cutoff      = input$enrichr_fdr_cutoff,
                                       db              = input$enrichr_db)
      
      rv$x2k_result <- run_x2k(rv$de_result, 
                            gene_direction = input$x2k_gene_direction, 
                            log2fc_cutoff  = input$x2k_abs_logFC_cutoff, 
                            fdr_cutoff     = input$x2k_fdr_cutoff,
                            dbs            = c('ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X','KEA_2015'))
    
  })
  
  output$exp_data_output <- renderDT(
    rv$exp,
    options = list(pageLength = 5, scrollX=TRUE)
  )
  
  output$colData_output <- renderDT(
    rv$colData,
    options = list(pageLength = 5, scrollX=TRUE)
  )
  
  output$pca_input <- renderUI({
    tagList(
      selectInput("pc_x", "Choose PC on x-axis", choices = 1:5, selected = 1),
      selectInput("pc_y", "Choose PC on y-axis", choices = 1:5, selected = 2),
      selectInput("pca_color", "Choose column to color", choices = unique(colnames(rv$colData))),
      actionButton("pca_button", "Go!")
    )
  })
  
  output$pca_scores_plot <- renderPlot({
    pca_plot(rv$pca_scores, input$pca_color, input$pc_x, input$pc_y)
  })
  
  output$de_result_table <- renderDT(
    rv$de_result,
    options = list(pageLength = 5, scrollX=TRUE)
  )
  
  output$de_volcano_plot <- renderPlot({
    volcano_plot(rv$de_result, 
                 log2fc_cutoff   = input$de_abs_logFC_cutoff, 
                 p_val_cutoff    = input$de_fdr_cutoff, 
                 ngenes_to_label = input$de_ngenes_to_label)
  })
  
  output$enrichr_dotplot <- renderPlot({
    plot_enriched_gene_sets(rv$egs, input$enrichr_ngs)
  })
  
  output$x2k_barplot <- renderPlot({
    plot_x2k_resuts(rv$x2k_result, input$x2k_nbars)
  })
}

shinyApp(ui, server)