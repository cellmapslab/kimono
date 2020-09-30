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

# Define UI for application
ui <- fluidPage(

    titlePanel("KiMONo"),

    sidebarPanel(
        fileInput("input_file", "Upload your different Omics Layers here", accept = c(".csv", ".tsv"), multiple= TRUE),
        fileInput("mapping_file", "Upload your Mapping Information here", accept = c(".csv", ".tsv"), multiple= TRUE),
        radioButtons("delimitor_input", "Delimiter",
                     choices = list("Space" = " ", "Tab" = "\t",
                                    "Comma" = ","),selected = " "),

        ),

    mainPanel(
        tableOutput("input_files"),
        tableOutput("mapping")
    )



)



# Define server logic required
server <- function(input, output) {

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
            x <- dim(read_delim(file$datapath[i], delim = input$delimitor_input2))
            x <- append(x, file$name[i], after = 0)
            print(x)
            input_tab <- rbind(input_tab,x )
        }
        rownames(input_tab) <- NULL
        colnames(input_tab) <- c("Filename", "Rows", "Columns")

        input_tab
    })
}

# Run the application
shinyApp(ui = ui, server = server)



