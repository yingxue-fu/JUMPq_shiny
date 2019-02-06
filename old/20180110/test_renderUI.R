## renderUI practice
library(shiny)

myDF <- data.frame(A = 1:4, B = 3:6, C = 6:9, D = 10:13)
ui = fluidPage(
  uiOutput("myList"),
  uiOutput("myNumbers")
)

server = function(input, output, session) {
  output$myList <- renderUI({
    selectInput("compName", "Company Name:", LETTERS[1:4])
  })
  output$myNumbers <- renderUI({
    selectInput("compNo", "Product Line:", myDF[, input$compName])
  })
}

shinyApp(ui, server)