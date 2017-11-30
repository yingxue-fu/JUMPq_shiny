library(shiny)

fluidPage (
  titlePanel('Quantification summary'),
  sidebarPanel(
    selectInput("level", label = "Select the quantification level",
                choices = c("Peptides", "Proteins"), 
                selected = "Peptides"),
    uiOutput("compSelection"),
    textInput("pval", label = "P-value threshold for differential expression (0 ~ 1)", value = 0.1),
    textInput("log2fc", label= "(+/- log2) Fold threshold for differential expression", value = 0.5),
    submitButton("Submit")
  ),

  ## Output visualization
  mainPanel(
    column(10, plotOutput("heatmap"))
  )
)















