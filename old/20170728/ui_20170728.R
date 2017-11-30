library(shiny)

fluidPage (
  titlePanel('Quantification summary'),
  sidebarPanel(
    selectInput("level", label = "Select the quantification level",
                choices = c("Unique peptides", "All peptides", "Unique proteins", "All proteins"), 
                selected = "Unique peptides"),
    textInput("fdr", label = "FDR threshold for differential expression (0 ~ 1)", value = 0.1),
    textInput("log2fc", label= "(+/- log2) Fold threshold for differential expression", value = 1),
    submitButton("Submit")
  ),

  ## Output visualization
  mainPanel(
    column(10, plotOutput("heatmap"))
  )
)















