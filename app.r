library(shiny)
library(Seurat)
library(ggplot2)

#Define Global Variables

load("data/day109.Rdata")

seurat.objects <- list(
  Day_20 = "Day_20",
  Day_109 = day109.seurat
)

feature.list <- list(
  Day_20 = "Day_20",
  Day_109 = day109.features
)

pal.list <- list(
  Day_20 = "Day_20",
  Day_109 = day109.palette
)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  #Suppress warnings
  tags$style(
    type="text/css",
    ".shiny-output-error { visibility: hidden; }",
    ".shiny-output-error:before { visibility: hidden; }"
  ),
  
  # App title ----
  titlePanel("Heatmap Customization"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel( width = 3,
      
      selectInput(
        "dataset",
        "Dataset",
        c("Day 109 Hair Cells" = "Day_109",
          "Day 20 Prosensory" = "Day_20"
          )
      ),
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "rmargin",
                  label = "Right Margin",
                  min = 0.05,
                  max = 0.5,
                  value = 0.15,
                  step = 0.01),
      
      #Adjust plot height
      sliderInput(inputId = "plotHeight",
                  label = "Plot Height",
                  min = 5,
                  max = 20,
                  value = 8.5,
                  step = 0.25),
      
      #Adjust plot width
      sliderInput(inputId = "plotWidth",
                  label = "Plot Width",
                  min = 5,
                  max = 20,
                  value = 11,
                  step = 0.25),
      
      sliderInput(inputId = "glabSize",
                  label = "Gene Label Size",
                  min = 0.2,
                  max = 6,
                  value = 2.0,
                  step = 0.1),
      
      sliderInput(inputId = "clabSize",
                  label = "Cluster Label Size",
                  min = 0.2,
                  max = 6,
                  value = 1.2,
                  step = 0.1),
      
      sliderInput(inputId = "legendScale",
                  label = "Legend Scale",
                  min = 0.05,
                  max = 1,
                  value = 0.5,
                  step = 0.05),
      
      checkboxInput("clusterLabels", "Cluster Labels", value = T),
      
      checkboxInput("legend", "Legend", value = T),
      
      downloadButton("HeatmapImage", "Download Plot as .png"),
      downloadButton("HeatmapEPS", "Download Plot as .eps file")
      
    ),
    
    
    
    # Main panel for displaying outputs ----
    mainPanel( width = 9,
      
      # Output: Histogram ----
      plotOutput(outputId = "Heatmap", inline = F)
      
    )
  )
)

server <- function(input, output) {
  
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  reactiveHeatmap <- reactive ({
    
    
    labeledHeat <- Seurat::DoHeatmap(seurat.objects[[input$dataset]],
                                     features = feature.list[[input$dataset]]$gene,
                                     group.colors = pal.list[[input$dataset]],
                                     label = input$clusterLabels,
                                     size = input$clabSize) +
      theme(plot.margin = unit(
                            c( 0.25,
                               input$rmargin,
                               0.25,
                               0.25), "in"),
            axis.text.y = element_text(size = input$glabSize),
            legend.key.size = unit(1.5 * input$legendScale, 'cm'), #change legend key size
            legend.key.height = unit(0.75  * input$legendScale, 'cm'), #change legend key height
            legend.key.width = unit(0.75 * input$legendScale, 'cm'), #change legend key width
            legend.title = element_text(size = 5.5 * input$legendScale), #change legend title font size
            legend.text = element_text(size = 5.5 * input$legendScale), #change legend text font size
            legend.spacing = unit(0.1 * input$legendScale, 'cm')
      )
    
    if (!input$legend) { labeledHeat <- labeledHeat + Seurat::NoLegend()}
    labeledHeat
  })
  
  
  
  reactiveHeight <- reactiveVal(500)
  observeEvent(input$plotHeight, {
    reactiveHeight( input$plotHeight*100 )
  })
  
  reactiveWidth <- reactiveVal(500)
  observeEvent(input$plotWidth, {
    reactiveWidth( input$plotWidth*100 )
  })
  
  
  observe({
    output$Heatmap <- renderPlot({
      
      reactiveHeatmap() },
      width = reactiveWidth(),
      height = reactiveHeight(),
      res = 300
    )
    
  })
  
  output$HeatmapImage <- downloadHandler(
    filename = function() { paste(input$dataset, "_Heatmap", '.png', sep='') },
    content = function(file) {
      ggsave(file,
             plot = reactiveHeatmap(),
             device = "png",
             units = "in",
             width = reactiveWidth()/300,
             height = reactiveHeight()/300,
             dpi = 300
      )
    }
  )
  
  output$HeatmapEPS <- downloadHandler(
    filename = function() { paste(input$dataset, "_Heatmap", '.eps', sep='') },
    content = function(file) {
      ggsave(file,
             plot = reactiveHeatmap(),
             device = "eps",
             units = "in",
             width = reactiveWidth()/300,
             height = reactiveHeight()/300,
             dpi = 300
      )
    }
  )
  
}

shinyApp(ui = ui, server = server)
