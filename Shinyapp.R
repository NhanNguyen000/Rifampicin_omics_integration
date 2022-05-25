library(shiny)
library(DT)
#library(tidyverse)
library(ggplot2)
load("movies.RData")


ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = "y", label = "Y-axis:",
                  choices = c("imdb_rating", "imdb_num_votes", "critics_score",
                              "audience_score", "runtime"),
                  selected = "audience_score"),
      selectInput(inputId = "x", label = "X-axis:",
                  choices = c("imdb_rating", "imdb_num_votes", "critics_score",
                              "audience_score", "runtime"),
                  selected = "critics_score"),
      sliderInput(inputId = "alpha", label = "Alpha",
                  min = 0, max = 1, value = 0.5),
      checkboxInput(inputId = "show_data", label = "Show data table", value = TRUE)
      ),
    mainPanel(plotOutput(outputId = "scatterplot"),
             plotOutput(outputId = "densityplot", height = 200),
              dataTableOutput(outputId = "moviestable"))
  )
)

server <- function(input, output, session) {
  output$scatterplot <- renderPlot({
    ggplot(data=movies, aes_string(x=input$x, y=input$y)) +
      geom_point(alpha = input$alpha)})
  output$densityplot <- renderPlot({
    ggplot(data = movies, aes_string(x=input$x)) +
      geom_density()})
  output$moviestable <- renderDataTable({
    if(input$show_data) {
      DT::datatable(data = movies %>% select(1:7),
                    options = list(pageLength=10), rownames = FALSE)}
  })
}
shinyApp(ui=ui, server = server)