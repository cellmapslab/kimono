#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(shiny)
library(data.table)
library(foreach)
library(sourcetools)
library(data.table)

#source("R/kimono.R")
#source("R/infer_sgl_model.R")
#source("R/utility_functions.R")

# Define UI for application
setwd("/Users/felicitastengen/Helmholtz_2/kimono/")
#source("20191121_kimono_test.R")

options(shiny.maxRequestSize = 30*1024^2)

ui <- pageWithSidebar(
    headerPanel("Kimono"),

    # Interface for the Sidebar
    sidebarPanel(

        # Sidebar for first Tab: Upload of Input Files
        conditionalPanel(condition="input.tabselected==1",
                         fileInput("input_file", "Upload your different Omics Layers here", accept = c(".csv", ".tsv"), multiple= TRUE),
                         radioButtons("delimitor_input", "Delimiter",
                                                         choices = list("Space" = " ", "Tab" = "\t","Comma" = ","),selected = " "),
                         uiOutput("selinput1")),

        # Sidebar for second Tab: Upload of Mapping Files
        conditionalPanel(condition="input.tabselected==2",
                         fileInput("mapping_file", "Upload your Mapping Information here", accept = c(".csv", ".tsv"), multiple= TRUE),
                                  radioButtons("delimitor_input_2", "Delimiter",
                                               choices = list("Space" = " ", "Tab" = "\t",
                                                              "Comma" = ","),selected = " "),),

        # Sidebar for third Tab: Choosing the Layers for each Mapping
        conditionalPanel(condition="input.tabselected==3", uiOutput("selmapping"),
                         uiOutput("selinput2"),
                         actionButton("do", "Add Mapping Info"),
                         actionButton("undo", "Delete current Mapping Info")),

        # Sidebar for fourth Tab: Main Settings for Kimono and Calculation
        conditionalPanel(condition="input.tabselected==4",
                         sliderInput("minfeatures", "Minimun Features:",
                                     min = 0, max = 20,
                                     value = 2),
                         actionButton("calculate", "Calculate Network")
                         )

    ),

    # Interface for the Main Page
    mainPanel(
        # Different Tabs that show different Outputs
        tabsetPanel(
            tabPanel("Input", value=1, tableOutput("input_files"), tableOutput("df_data_out")),
            tabPanel("Mapping", value=2,tableOutput("mapping_files") ),
            tabPanel("Metadata", value=3, tableOutput("user_data")),
            tabPanel("Calculation", value = 4),

            id = "tabselected"
        )
    )
)





server <- function(input, output) {


    # observes the Input Files that are uploaded by the User
    # Each file is read into a Dataframe and saved in the reactive Values Object
    layers <- reactiveValues()

    observeEvent( input$input_file ,{

        filenames <- input$input_file$name
        for(i in 1:length(filenames)){
            layers[[filenames[i]]] <- as.data.table(read_delim(input$input_file$datapath[i], delim = input$delimitor_input))
            layers[[filenames[i]]] <- layers[[filenames[i]]][,-1,with=FALSE]
        }

        print(isolate(reactiveValuesToList(layers)))


    })


    # observes the Input Files that are uploaded by the User
    # Each file is read into a Dataframe and saved in the reactive Values Object
    mappings <- reactiveValues()
    observeEvent(input$mapping_file,{

        filenames <- input$mapping_file$name
        for(i in 1:length(filenames)){
            mappings[[filenames[i]]] <- as.data.table(read_delim(input$mapping_file$datapath[i], delim = input$delimitor_input_2))
        }

        print(isolate(reactiveValuesToList(mappings)))

    })








    # Iterates over the input Files in the reactive layers object and calculates a new dataframe
    # For each uploaded File it shows the filename and the dimensions when read ito a Dataframe
    # this table is saved in the reactive Value whichfiles object and also directly generated as output
    # by the renderTable() Function
    whichfiles <- reactiveValues()
    output$input_files <- renderTable({
         input_tab <- NULL

         if(length(names(layers) >0)){
         for (i in 1:length(names(layers))){

             x <- dim(layers[[names(layers)[i]]])
             x <- append(x, names(layers)[i], after = 0)
             print(x)
             input_tab <- rbind(input_tab,x )

         }
         rownames(input_tab) <- NULL
         colnames(input_tab) <- c("Filename", "Rows", "Columns")
         whichfiles$input <- input_tab
         input_tab
         }

     })



    # Iterates over the Mapping Files in the reactive layers object and calculates a new dataframe
    # For each uploaded File it shows the filename and the dimensions when read ito a Dataframe
    # this table is saved in the reactive Value whichfiles object and also directly generated as output
    # by the renderTable() Function

     output$mapping_files <- renderTable({
         mapping_tab <- NULL

         if(length(names(mappings) >0)){
             for (i in 1:length(names(mappings))){

                 x <- dim(mappings[[names(mappings)[i]]])
                 x <- append(x, names(mappings)[i], after = 0)
                 print(x)
                 mapping_tab <- rbind(mapping_tab,x )

             }
             rownames(mapping_tab) <- NULL
             colnames(mapping_tab) <- c("Filename", "Rows", "Columns")
             whichfiles$mapping <- mapping_tab
             mapping_tab
         }

     })

    #' @TODO BUG if only one file is uploaded than the dropdown does not work - Exception should be introduced

    # Dropdown Selection of the Input File Names that were uploaded by the user that can be accesed in UI
     output$selinput1 <- renderUI({
         selectInput(inputId = "selection1",
                     label = "Select the Main Layer",
                     selected = whichfiles$input[,1],
                     choices = whichfiles$input[,1])
     })

     # Dropdown Selection of the Input File Names that were uploaded by the user that can be accesed in UI
     output$selinput2 <- renderUI({
         selectInput(inputId = "selection2",
                     label = "Select Layer",
                     selected = whichfiles$input[,1],
                     choices = whichfiles$input[,1])
     })

     # Dropdown Selection of the Mapping File Names that were uploaded by the user that can be accesed in UI
     output$selmapping <- renderUI({
         selectInput(inputId = "mapping",
                     label = "Select Mapping File",
                     selected = whichfiles$mapping[,1],
                     choices = whichfiles$mapping[,1])
     })




     # observes the Button and adds the information from the dropdown menus to a dataframe that is saved in the
     # meta reactive Value object
     meta <- reactiveValues(info = NULL)

     observeEvent(input$do, {

         meta$info <- rbind(meta$info, c( input$mapping, input$selection1, input$selection2) )
         print(meta$info)
         colnames(meta$info) <- c("mapping_file", "Layer_1", "Layer_2")
         output$user_data <- renderTable(meta$info)

     })

     # Clears the whole DataFrame
     observeEvent(input$undo,{
         meta$info <- NULL

     })

     #' @TODO Input to Kimono and upload the resulting DataFrame to the User Interface
     observeEvent(input$calculate,{

         mymeta <- as_tibble(isolate(meta$info))

        print(names(layers))
        #print(mymeta$Layer_2)
        mymeta$layernumber <- match(mymeta$Layer_2, names(layers))
        mymeta <- mymeta[match(names(mappings), mymeta$mapping_file),]

        metainfo <- data.frame('ID' = names(mappings), 'main_to' = mymeta$layernumber)
        omicdf <- isolate(reactiveValuesToList(layers))
        mappingdf <- isolate(reactiveValuesToList(mappings))
        mininput <- input$minfeatures
        mainlayer <- match(input$selection1, names(layers))

        print(metainfo)
        print(omicdf)
        print(mappingdf)
        print(mainlayer)
        print(mininput)


        print(kimono(input_list = omicdf, mapping_list = mappingdf,metainfo = metainfo, main_layer = main_layer, min_features = mininput ))
     })

 }

# Run the application
shinyApp(ui = ui, server = server)



