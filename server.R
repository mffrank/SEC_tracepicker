library(ggplot2)
library(data.table)
library(reshape)
library(plotly)
library(DT)

# library(Cairo)   # For nicer ggplot2 output when deployed on Linux
# source("plot.complexFeatures_peptides.R")
source("tracesMethods.R")

# Load data
pepTraces.annotated.consec3.SibPepCorr.filtered <- readRDS("data/pepTracses.annotated.consec3.SibPepCorr.filtered.all.rda")
features <- readRDS("data/manual_annotation_all_revised.rds")
hyp_names <- paste(features$protein_id,features$apex,sep = "_")
features$id <- hyp_names

instance <- round(runif(1)*10^7,0)

server <- function(input, output, session) {
  annotations <- reactiveValues()
  annotations$dt <- data.table(FeatureName = numeric(0),
                               leftBoundary = numeric(0),
                               rightBoundary = numeric(0),
                               Confidence = numeric(0),
                               Apex = numeric(0),
                               Subunits = numeric(0))
  ranges <- reactiveValues(x = NULL, y = NULL)
  subunits <- reactiveValues(sub = NULL)
  # feature <- reactive(features[id == input$select])
  
  ## Observe the Button Events ------------------
  ###################################################
  
  plot.data <- eventReactive(input$select,{
    feature <- features[id == input$select]
    ranges$x <- c(feature$left_pp,feature$right_pp)

    peptides <- pepTraces.annotated.consec3.SibPepCorr.filtered$trace_annotation$id[pepTraces.annotated.consec3.SibPepCorr.filtered$trace_annotation$protein_id == feature$protein_id]
    pepTraces_sel <- subset(pepTraces.annotated.consec3.SibPepCorr.filtered,trace_ids=peptides)
    toLongFormat(pepTraces_sel$traces)
  })
  
  observeEvent(input$prev, {
    pos <- which(input$select == hyp_names)
    if(pos != 1){
      updateSelectInput(session, "select", selected = hyp_names[pos-1])
      subunits$sub = NULL
     }
  })
  
  observeEvent(input$forw, {
    pos <- which(input$select == hyp_names)
    if(pos != length(hyp_names)){
      updateSelectInput(session, "select", selected = hyp_names[pos+1])
      subunits$sub = NULL
    }
  })
  
  observeEvent(input$save_sure, {
    if(!is.null(subunits$sub)){
      newLine <- cbind(features[id == input$select,.(FeatureName = protein_id,
                                                     Apex = apex,
                                                     leftBoundary = left_pp,
                                                     rightBoundary = right_pp,
                                                     Confidence = confidence)],
                       Subunits = paste(subunits$sub, collapse = ";"))
      isolate(annotations$dt <- rbind(newLine,annotations$dt))
      subunits$sub <- NULL
    }
    })
  
  observeEvent(input$unzoom, {
    ranges$y <- NULL
    # ranges$y <- NULL
  })
  
  observeEvent(input$saveall, {
    saveRDS(isolate(annotations$dt), file = paste0("data/Manual_annotation",gsub("-","_",Sys.Date()),"_",instance,".rda"))
  })
  
  observeEvent(input$delete, {
    if(!is.null(input$view_rows_selected)){
      annotations$dt<- annotations$dt[-input$view_rows_selected]
    }
  })
  
  observeEvent(input$keypress[1],{
    if(input$keypress[1] == "3" | input$keypress[1] == "ArrowRight"){
      pos <- which(input$select == hyp_names)
      if(pos != length(hyp_names)){
        updateSelectInput(session, "select", selected = hyp_names[pos+1])
        subunits$sub = NULL
      } 
    } else if(input$keypress[1] == "2" | input$keypress[1] == "ArrowLeft"){
      pos <- which(input$select == hyp_names)
      if(pos != 1){
        updateSelectInput(session, "select", selected = hyp_names[pos-1])
        subunits$sub = NULL
      }
    } else if(input$keypress[1] == "1"){
      if(!is.null(subunits$sub)){
        newLine <- cbind(features[id == input$select,.(FeatureName = protein_id,
                                                       Apex = apex,
                                                       leftBoundary = left_pp,
                                                       rightBoundary = right_pp,
                                                       Confidence = confidence)],
                         Subunits = paste(subunits$sub, collapse = ";"))
        isolate(annotations$dt <- rbind(newLine,annotations$dt))
        subunits$sub <- NULL
      }
    } 
    # print(input$keypress)
  })
  
  observeEvent(event_data("plotly_click"),{
    sel <- event_data("plotly_click")$key
    if(sel %in% subunits$sub){
      subunits$sub <- subunits$sub[!(subunits$sub == sel)]
    } else {
      subunits$sub <- c(subunits$sub, sel)#gsub("^;","",paste(subunits$sub, sel$key, sep = ";"))
    }
    print(subunits$sub)
  })
  
  observeEvent(event_data("plotly_selected"),{
    sel <- unique(event_data("plotly_selected")$key)
    subunits$sub <- c(subunits$sub, setdiff(sel,subunits$sub))
  })
  ## Observe the Mouse Events in Plots ---------------------------
  ###############################################################
  
  observeEvent(input$plot2_dblclick, {
    brush <- input$plot2_brush
    if (is.null(brush)) {
      ranges$y <- NULL
    }
  })

  observeEvent(input$plot2_brush, {
    brush <- input$plot2_brush
    ranges$y <- c(round(brush$ymin,0), round(brush$ymax,0))
  })
  
  
  ## Output the Plots and Data table -------------------------
  ###############################################
  
  
  output$plot1 <- renderPlotly({
    # print(apex_sel())
    p <- plot.traces(plot.data(),ranges,apex_man = list(sel= features[id == input$select]$apex),
                     plot=FALSE, ledgend = FALSE, log =(input$log == "Log"))
    ggplotly(p) %>% layout(dragmode= "select", yaxis = list(title = "log(intensity)"))
    
  })
  
  output$plot2 <- renderPlot({

    p <- plot.traces(plot.data(),ranges = list(x = ranges$x, y = NULL), plot=TRUE,
                            ledgend = FALSE,log = FALSE)
    p
  })
  
  output$view <-DT::renderDataTable(annotations$dt)
  
  # ## Used for debugging, no longer needed-------------------
  output$click_info <- renderPrint({
    # Because it's a ggplot2, we don't need to supply xvar or yvar; if this
    # were a base graphics plot, we'd need those.
    # nearPoints(pepTraces.annotated.consec3.SibPepCorr.filtered, input$plot1_click, addDist = TRUE)
    # s <- event_data("plotly_click")$key
    # if (length(s) == 0) {
    #   "Click on a trace in the plot to select it"
    # } else {
    #   cat("You selected: \n\n")
    #   as.list(s)
    # }
    subunits$sub
  })
  
  # output$brush_info <- renderPrint({
  #   s <- event_data("plotly_selected")
  #   if (length(s) == 0) {
  #     "Click on a cell in the plot to display values"
  #   } else {
  #     cat("You selected: \n\n")
  #     as.list(s)
  #   }
  #   })
  #----------------------------------
  
}
