## For information, https://github.com/rstudio/d3heatmap

rm(list = ls())
# library(d3heatmap)
# library(shiny)
# ui <- fluidPage(
#   h1("A heatmap demo"),
#   sliderInput("pixels", "size", value = 400, min = 100, max = 1000),
#   selectInput("palette", "Palette", c("YlOrRd", "RdYlBu", "Greens", "Blues")),
#   checkboxInput("cluster", "Apply clustering"),
#   uiOutput("dynamic")
# )
# 
# server <- function(input, output) {
#   output$heatmap <- renderD3heatmap({
#     d3heatmap(
#       scale(mtcars),
#       colors = input$palette,
#       dendrogram = if (input$cluster) "both" else "none"
#     ) })
#   
#   output$dynamic <- renderUI({
#     
#     d3heatmapOutput("heatmap", height = paste0(input$pixels, "px"))
#   })
# }
# 
# shinyApp(ui, server)


library(d3heatmap)
library(shiny)

ui <- fluidPage(
  h1("A heatmap demo"),
  selectInput("palette", "Palette", c("YlOrRd", "RdYlBu", "Greens", "Blues")),
  checkboxInput("cluster", "Apply clustering"),
  d3heatmapOutput("heatmap")
)

server <- function(input, output, session) {
  output$heatmap <- renderD3heatmap({
    d3heatmap(
      scale(mtcars),
      colors = input$palette,
      dendrogram = if (input$cluster) "both" else "none"
    )
  })
}

shinyApp(ui, server)
