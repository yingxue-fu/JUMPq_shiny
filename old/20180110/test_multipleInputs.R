library(shiny)

ui = fluidPage(
  
  headerPanel("Example"),
  
  sidebarPanel(
    div(class="row-fluid",
        div(class="span3",
            checkboxGroupInput("one", "One:",letters[1:26])),
        div(class="span3",
            checkboxGroupInput("two", "Two:",letters[1:26]))),
    
    selectInput("three", label="Three", choices=letters, multiple=TRUE)
    
  ),
  mainPanel()
)

server = function(input, output){}

shinyApp(ui = ui, server = server)