rm(list = ls())

library(shiny)

ui = fluidPage(
  sliderInput(inputId = "num",
              label = "Choose a number",
              value = 25, min = 1, max = 100),
  plotOutput("hist")
)

server = function(input, output) {

  ## This causes an error of reactivity since the reactive input is changed 
  ## outside render function
  bin = input$num 
  output$hist = renderPlot({
    hist(rnorm(bin))
  })
    
  ## All the reactivity should be implemented within "render~" function
  output$hist = renderPlot({
    bin = input$num 
    hist(rnorm(bin))
  })
}

shinyApp(ui = ui, server = server)