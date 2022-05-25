library(shiny)
library(visNetwork)
load("./outcome/visNetwork.RData")

ui <- fluidPage(
  titlePanel("Uploading Files"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose CSV file", multiple = FALSE,
                accept = c("test/csv", "text/comma-separated-values,text/plain", ".csv")),
      tags$hr(),
      checkboxInput("header", "Header", TRUE),
      radioButtons("sep", "Separator",
                   choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
                   selected = ","),
      radioButtons("quote", "Quote",
                   choices = c(None = "", "Double Quote" = '"', "Single Quote" = "'"),
                   selected = '"'),
      radioButtons("disp", "Display",
                   choices = c(Head = "head", All = "all"),
                   selected = "head")
    ),
  mainPanel(
    tabsetPanel(
      tabPanel("Table", tableOutput("contents")),
      tabPanel("Network", visNetworkOutput("network_proxy_nodes", height = "600px"))
      )
    )
  )
)

server <- function(input, output, session) {
  output$contents <- renderTable({
    req(input$file1)
    tryCatch({
      df <- read.csv(input$file1$datapath, header = input$header,
                     sep = input$sep, quote = input$quote)
    }, 
    error = function(e) {stop(safeError(e))})
    if(input$disp == "head") {return(head(df))} else {return(df)}
  })
  
  output$network_proxy_nodes <- renderVisNetwork({
    visnet <- visNetwork(nodes.vis, edges.vis) %>% visGroups(groupname = "TF", color = "orange") %>%
      visGroups(groupname = "target", color = "lightblue") %>%
      visLegend() %>% visInteraction(navigationButtons = TRUE)
    visOptions(visnet, 
               nodesIdSelection = TRUE, highlightNearest = TRUE,
               selectedBy = "group")
  })
}

shinyApp(ui, server)