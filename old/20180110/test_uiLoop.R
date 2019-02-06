library(shiny)

server = function(input,output) {
  observe({
    if (input$test == 0) 
      return()
    isolate({
      output$value <-renderTable({
        matrix_input <- list()
        for(j in 1:(input$obs+1)){matrix_input[j] <- paste0("<input id='element",j,"_", 1:1,
                                                            "' class='shiny-bound-input span7' type='number' value=''>")}
        matrix <- data.frame(matrix_input) 
        colnames(matrix)<- c("Size n",sapply(1:input$obs,function(z)paste("Cat",z)))
        matrix}, sanitize.text.function = function(x) x)
    })
  })
  output$summary <- renderPrint({summary("'element$value")})
}


ui = fluidPage (
  headerPanel(title ="Elicitación de Expertos"),
  sidebarPanel(
    tags$head(
      tags$style(type="text/css", "select { width: 100px; }"),
      tags$style(type="text/css", "textarea { max-width: 185px; }"),
      tags$style(type="text/css", ".jslider { max-width: 200px; }"),
      tags$style(type='text/css', ".well { max-width: 470px; }"),
      tags$style(type='text/css', ".span4 { max-width: 830px; }")
    ), 
    numericInput("obs", "Número de Categorías:", 2, min =2, max =10),
    h3(actionButton("test","Crear")),
    
    conditionalPanel(condition = "input.density == true",sliderInput(inputId = "bw_adjust", label = "Bandwidth Adjustment",
                                                                     min = 0.2, max = 2, value = 1, step = 0.2))
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Plot", tableOutput("value")),
      tabPanel("Summary",verbatimTextOutput(outputId = "summary"))
    ))
)

shinyApp(ui, server)