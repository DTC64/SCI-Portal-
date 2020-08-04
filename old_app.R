library('shiny')
library('ggplot2')
library('cowplot')
library('dplyr')
library('ggrepel')

source('helpers.R')

LR_results <- readRDS('data/LR_results.rds')
LR_results$Pair_name <- factor(LR_results$Pair_name, levels = sort(levels(LR_results$Pair_name)))
LR_results$Ligand <- factor(LR_results$Ligand, levels = sort(levels(LR_results$Ligand)))
LR_results$Receptor <- factor(LR_results$Receptor, levels = sort(levels(LR_results$Receptor)))

sci_umap <- readRDS('data/sci_umap.rds')

ui <- fluidPage(
  tabsetPanel(
    tabPanel(
      title = 'Ligand-Receptor Interactions after SCI',
      fluid = TRUE,
      h3('LR Plotting options'),
      fluidRow(
        column(
          width = 3, 
          selectInput(
            inputId = 'Pair_name',
            label = 'Select specific Ligand-Receptor Pairs:',
            choices = {l <- as.list(levels(LR_results$Pair_name)); names(l) <- levels(LR_results$Pair_name); l},
            multiple = TRUE,
            selectize = TRUE
          ),
          selectInput(
            inputId = 'Ligand',
            label = 'Include all pairs with the following Ligands (leave blank for all ligands):',
            choices = {l <- as.list(levels(LR_results$Ligand)); names(l) <- levels(LR_results$Ligand); l},
            multiple = TRUE,
            selectize = TRUE
          ),
          selectInput(
            inputId = 'Receptor',
            label = 'Include all pairs with the following Receptors (leave blank for all receptors):',
            choices =  {r <- as.list(levels(LR_results$Receptor)); names(r) <- levels(LR_results$Receptor); r},
            multiple = TRUE,
            selectize = TRUE
          )
        ),
        column(
          width = 3,
          selectInput(
            inputId = 'Ligand_Cell',
            label = 'Search ligand expression in these cells (leave blank for all cells):',
            choices =  {lc <- as.list(c(names(cell), levels(LR_results$Ligand_Cell))); names(lc) <- c(names(cell), levels(LR_results$Ligand_Cell)); lc},
            multiple = TRUE,
            selectize = TRUE
          ),
          selectInput(
            inputId = 'Receptor_Cell',
            label = 'Search receptor expression in these cells (leave blank for all cells):',
            choices = {rc <- as.list(c(names(cell), levels(LR_results$Receptor_Cell))); names(rc) <- c(names(cell),levels(LR_results$Receptor_Cell)); rc},
            multiple = TRUE,
            selectize = TRUE
          )
        ),
        column(
          width = 3,
          checkboxGroupInput(
            inputId = 'Time',
            label = 'Time-points:',
            choices = {t <- as.list(levels(LR_results$Time)); names(t) <- levels(LR_results$Time); t},
            selected = levels(LR_results$Time)
          ),
          h5(strong('Interaction Enrichment')),
          checkboxInput(
            inputId = 'significant.only',
            label = 'Display P-val < 0.05 only?',
            value = TRUE
          )
        ),
        column(
          width = 3,
          radioButtons(
            inputId = 'organize',
            label = 'Group Y-axis by:',
            choices = list('Ligand Cell' = 'Ligand_Cell',
                           'Receptor Cell' = 'Receptor_Cell'),
            selected = 'Ligand_Cell'
          ),
          radioButtons(
            inputId = 'sort.order',
            label = 'Sort LR pairs by:',
            choices = list('Ligand' = 'Ligand',
                           'Receptor' = 'Receptor'),
            selected = 'Ligand'
          ),
          actionButton(
            inputId = 'ApplyLRChange', 
            label = 'Apply changes'
          )
        )
      ),
      uiOutput('LRplot.ui')
    )
  ),
  tabsetPanel(
    tabPanel(
      title = 'UMAP after SCI',
      fluid = TRUE,
      h3('UMAP Plotting options'),
      sidebarLayout(
        sidebarPanel(
          width = 4,
          selectInput(
            inputId = 'color.by', 
            label = 'Color cells by (gene names give expression):',
            choices = {x <- as.list(colnames(sci_umap)); names(x) <- colnames(sci_umap); x}
          ),
          checkboxInput(
            inputId = 'split.time',
            label = 'Split data by timepoint?',
            value = TRUE
          ),
          br(),
          h4('Reference UMAP'),
          img(src = 'Fig1a_Heterogeneity_After_SCI_1.png', width = 450, height = 400),
          br(),
          h4('All sub-clusters'),
          img(src = 'SCI_all_clusters.png', width = 450, height = 400)
        ),
        mainPanel(
          uiOutput('UMAPplot.ui')
        )
      )
    )
  )
)


server <- function(input, output) {
  
  tmp <- eventReactive(
    eventExpr = input$ApplyLRChange,
    valueExpr = {
      selectLR(
        score.results.df = LR_results,
        timepoints = input$Time,
        receptor.cell = input$Receptor_Cell,
        ligand.cell = input$Ligand_Cell,
        ligand = input$Ligand,
        receptor = input$Receptor,
        pair.name = input$Pair_name,
        significant.only = input$significant.only,
        organize = input$organize,
        sort.order = input$sort.order
      )
    }
  )
  
  h <- eventReactive(
    eventExpr = input$ApplyLRChange,
    valueExpr = {
      1.5*length(unique(tmp()$data[[tmp()$y.axis]])+3) * 2*(length(unique(tmp()$data[[tmp()$y.facet]]))+10) * (length(unique(tmp()$data$Time))+3)
    }
  )
  
  output$LRplot <- renderPlot({
    obj <- tmp()
    do.call(
      what = obj$style$p,
      args = list(
        'score.results.df' = obj$data, 
        'x.axis' = obj$x.axis,
        'y.axis' = obj$y.axis,
        'y.facet' = obj$y.facet,
        'switch' = obj$switch,
        'dims' = obj$style$dims
      )
    )
  })
  
  output$LRplot.ui <- renderUI({
    plotOutput('LRplot', height = h())
  })
  
  umap <- reactive({
    selectUMAP(
      df = sci_umap,
      color.by = input$color.by
    )
  })
  
  output$UMAPplot <- renderPlot({
    obj <- umap()
    UMAP_LR(obj, color.by = input$color.by, split.time = input$split.time)
  })
  
  output$UMAPplot.ui <- renderUI({
    plotOutput('UMAPplot', height = 1050)
  })
  
}

shinyApp(ui, server)

# library(shiny)
# 
# ui <- fluidPage(
#   titlePanel('My Shiny App'),
#   sidebarLayout(
#     
#     sidebarPanel(
#       h3('Installation'),
#       p('Shiny is available on blah blah from you console:', code('www.google.com')),
#       br(),
#       br(),
#       img(src = 'download.jfif', height = 140, width = 140, align = 'center'),
#       br(),
#       p('Miami Project', align = 'center')
#     ),
#     
#     mainPanel(
#       h1('Introducing Shiny'),
#       p('Shiny is a new application blah blah it is totally', em('amazinggggg')),
#       br(),
#       p('more examples can be found here:', a('www.google.com'))
#     )
#   )
# )
# 
# server <- function(input, output, session) {
#   
# }
# 
# shinyApp(ui, server)
