library(shiny)
library(ggplot2)
library(dplyr)
library(DT)

load("movies.RData")

server <- function(input, output, session) {
  output$moviestable <- renderDataTable({
    req(input$n)
    movies_sample <- movies %>% sample_n(input$n) %>% select(title:studio)
    DT::datatable(data = movies_sample, options = list(pageLength =10), rownames = FALSE)
  })
}

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      HTML(paste("Enter a value between 1 and", "651")),
      numericInput(inputId = "n", label = "Sample size:",
                   value = 3, step = 10)
    ),
    mainPanel(DT::dataTableOutput(outputId = "moviestable"))
  )
)

shinyApp(ui = ui, server = server)