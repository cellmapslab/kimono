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

# Define UI for application

options(shiny.maxRequestSize = 30*1024^2)

ui <- pageWithSidebar(
    headerPanel("Kimono"),
    sidebarPanel(

        ## conditionalPanel() functions for selected tab
        conditionalPanel(condition="input.tabselected==1",
                         fileInput("input_file", "Upload your different Omics Layers here", accept = c(".csv", ".tsv"), multiple= TRUE),
                         radioButtons("delimitor_input", "Delimiter",
                                                         choices = list("Space" = " ", "Tab" = "\t","Comma" = ","),selected = " ")),


        conditionalPanel(condition="input.tabselected==2",
                         fileInput("mapping_file", "Upload your Mapping Information here", accept = c(".csv", ".tsv"), multiple= TRUE),
                                  radioButtons("delimitor_input_2", "Delimiter",
                                               choices = list("Space" = " ", "Tab" = "\t",
                                                              "Comma" = ","),selected = " "),),


        conditionalPanel(condition="input.tabselected==3", uiOutput("selmapping"),
                         uiOutput("selinput1"),
                         uiOutput("selinput2"),
                         actionButton("do", "Add Mapping Info"))

    ),
    mainPanel(
        # recommend review the syntax for tabsetPanel() & tabPanel() for better understanding
        # id argument is important in the tabsetPanel()
        # value argument is important in the tabPanle()
        tabsetPanel(
            tabPanel("Input", value=1, tableOutput("input_files"), tableOutput("df_data_out")),
            tabPanel("Mapping", value=2,tableOutput("mapping") ),
            tabPanel("Metadata", value=3, tableOutput("user_data")),
            id = "tabselected"
        )
    )
)





server <- function(input, output) {

    mynames <- c("a", "b")
    values <- reactiveValues(df_data = NULL)
    observeEvent(input$input_file, {



        values[[mynames[1]]] <- read_delim(input$input_file$datapath[1], delim = " ")
    })

    output$df_data_out <- renderTable(values[[mynames[1]]])




    whichfiles <- reactiveValues()

     output$input_files <- renderTable({
         file <- input$input_file
         print(file$name)

         ext <- tools::file_ext(file$datapath)


         input_tab <- NULL
         req(file)
         validate(need(ext == c("csv","tsv"), "Please upload a csv or tsv file"))


         for (i in 1:length(file$name)){
             x <- dim(read_delim(file$datapath[i], delim = input$delimitor_input))
             x <- append(x, file$name[i], after = 0)
             print(x)
             input_tab <- rbind(input_tab,x )
         }
         rownames(input_tab) <- NULL
         colnames(input_tab) <- c("Filename", "Rows", "Columns")

        whichfiles$input <- input_tab
        input_tab
     })



     output$mapping <- renderTable({
         file <- input$mapping_file
         print(file$name)

         ext <- tools::file_ext(file$datapath)


         input_tab <- NULL
         req(file)
         validate(need(ext == c("csv","tsv"), "Please upload a csv or tsv file"))


         for (i in 1:length(file$name)){
             x <- dim(read_delim(file$datapath[i], delim = input$delimitor_input_2))
             x <- append(x, file$name[i], after = 0)
             print(x)
            input_tab <- rbind(input_tab,x )
         }
         rownames(input_tab) <- NULL
         colnames(input_tab) <- c("Filename", "Rows", "Columns")
         whichfiles$mapping <- input_tab
         input_tab
     })




     output$selinput1 <- renderUI({
         selectInput(inputId = "selection1",
                     label = "Select Layer",
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


     observeEvent(input$do, {
         user_data <- data.frame(
             "mapping_file" =input$mapping,
             "Layer_1" =input$selection1,
             "Layer_2" =input$selection2,
             stringsAsFactors = FALSE
         )

         # ptdata = merge(pub_data,user_data,by = c("gender","agemo","racethn"), all = TRUE)
         # names(ptdata)<-tolower(names(ptdata))
         # ptdata<-ptdata %>%
         #     mutate(corrected_nutrient = nutrient/coef)

         output$user_data <- renderTable(user_data)
         #output$FinalData <- DT::renderDataTable(DT::datatable(ptdata))
     })


 }

# Run the application
shinyApp(ui = ui, server = server)



