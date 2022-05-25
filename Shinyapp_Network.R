library(shiny)
library(visNetwork)
load("./outcome/visNetwork.RData")

ui <- fluidPage(
  fluidRow(
    column(width = 2,  tableOutput("contents")),
    column(width = 8, visNetworkOutput("network_proxy_nodes", height = "600px"))
  )
)

server <- function(input, output, session) {

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