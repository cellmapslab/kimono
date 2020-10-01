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
library(kimono)

# Define UI for application

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
            tabPanel("Mapping", value=2,tableOutput("mapping") ),
            tabPanel("Metadata", value=3, tableOutput("user_data")),
            tabPanel("Calculation", value = 4),

            id = "tabselected"
        )
    )
)





server <- function(input, output) {



    layers <- reactiveValues()

    toListen <- reactive({
        list(input$input_file,input$delimitor_input)
    })


    observeEvent( input$input_file ,{

        filenames <- input$input_file$name
        for(i in 1:length(filenames)){
            layers[[filenames[i]]] <- read_delim(input$input_file$datapath[i], delim = input$delimitor_input)
        }

        print(isolate(reactiveValuesToList(layers)))


    })

    mappings <- reactiveValues()
    observeEvent(input$mapping_file,{

        filenames <- input$mapping_file$name
        for(i in 1:length(filenames)){
            mappings[[filenames[i]]] <- read_delim(input$mapping_file$datapath[i], delim = input$delimitor_input_2)
        }

        print(isolate(reactiveValuesToList(mappings)))

    })







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




     output$mapping <- renderTable({


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
         # file <- input$mapping_file
         # print(file$name)
         #
         # ext <- tools::file_ext(file$datapath)
         #
         #
         # input_tab <- NULL
         # req(file)
         # validate(need(ext == c("csv","tsv"), "Please upload a csv or tsv file"))
         #
         #
         # for (i in 1:length(file$name)){
         #     x <- dim(read_delim(file$datapath[i], delim = input$delimitor_input_2))
         #     x <- append(x, file$name[i], after = 0)
         #     print(x)
         #    input_tab <- rbind(input_tab,x )
         # }
         # rownames(input_tab) <- NULL
         # colnames(input_tab) <- c("Filename", "Rows", "Columns")
         # whichfiles$mapping <- input_tab
         # input_tab
     })




     output$selinput1 <- renderUI({
         selectInput(inputId = "selection1",
                     label = "Select the Main Layer",
                     selected = whichfiles$input[,1],
                     choices = whichfiles$input[,1])
     })

     output$selinput2 <- renderUI({
         selectInput(inputId = "selection2",
                     label = "Select Layer",
                     selected = whichfiles$input[,1],
                     choices = whichfiles$input[,1])
     })


     output$selmapping <- renderUI({
         selectInput(inputId = "mapping",
                     label = "Select Mapping File",
                     selected = whichfiles$mapping[,1],
                     choices = whichfiles$mapping[,1])
     })




     meta <- reactiveValues(info = NULL)

     observeEvent(input$do, {

         meta$info <- rbind(meta$info, c( input$mapping, input$selection1, input$selection2) )
         print(meta$info)
         colnames(meta$info) <- c("mapping_file", "Layer_1", "Layer_2")
         output$user_data <- renderTable(meta$info)

     })

     observeEvent(input$undo,{
         meta$info <- NULL

     })

     observeEvent(input$calculate,{






     })

 }

# Run the application
shinyApp(ui = ui, server = server)



