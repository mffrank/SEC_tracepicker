
library(plotly)

jscode <- '
$(function() {
$(document).keypress(function(e) {
Shiny.onInputChange("keypress", [e.key, Math.random()]);

});
});
'

# Load data
features <- readRDS("data/manual_annotation_all_revised.rds")
hyp_names <- paste(features$protein_id,features$apex,sep = "_")

ui <- fluidPage(
  tags$head(tags$script(HTML(jscode))),
  titlePanel("SEC-WATH MS Tracepicker"),
  fluidRow(
    column(4,
      wellPanel(
        h3("Select Feature"),
        br(),
        selectizeInput("select", label = h5("UniprotID"), 
                    choices = hyp_names,
                    options = list(maxOptions=10000)),
        actionButton("prev","",icon = icon("arrow-left")),
        actionButton("forw","",icon = icon("arrow-right")),
        br(),
        br(),
        fluidRow(
          # actionButton("zoom",label = "Zoom"),
          actionButton("unzoom",label = "Reset Zoom")
        ),
        br(),
        fluidRow(
          actionButton("save_sure", label = "Save Traces",icon = icon("thumbs-o-up"))
          # ,
          # actionButton("save", label = "Save Feature",icon = icon("thumbs-o-down"))
        ),
        actionButton("saveall", label = "Save Data Frame")
      ),
      plotOutput("plot2",
                 dblclick = "plot2_dblclick",
                 brush = brushOpts(id = "plot2_brush",
                                   direction = "y",
                                   resetOnNew = TRUE),
                 height = 200
                 ),
      h4("How to annotate Features"),
      tags$ul(
        tags$li("Use small plot to zoom to desired region"), 
        tags$li("Drag your mouse over the Traces to select multiple ones"), 
        tags$li("Click to select single traces"),
        tags$li("Click on selected trace to deselect it"),
        tags$li("Save the feature by pressing 'Save Traces'"),
        tags$li("Delete wrong entries by selecting them and pressing delete"),
        tags$li("Press Save Data Frame often, to save a copy of the table")
      )
      
    ),
    column(8,
      # conditionalPanel(condition = "$('html').hasClass('shiny-busy')",
      #       h4("The data is loading...")
      # ),
      radioButtons("log",label = NULL, choices = c("Log","Continous"),inline = TRUE),
      plotlyOutput("plot1"),
      
      # ## Debugging
      fluidRow(
        column(width = 12,
               verbatimTextOutput("click_info")
        )
      # ,
      # column(width = 6,
      #        verbatimTextOutput("brush_info")
      # )
      ),
    fluidRow(
        div(actionButton("delete", "Delete",icon = icon("remove-sign",lib="glyphicon")), style="float:right")
      ),
      # div(style="display:inline-block",actionButton("delete", "Delete"), style="float:right"),
      br(),
      DT::dataTableOutput("view")
    )
  )
)
