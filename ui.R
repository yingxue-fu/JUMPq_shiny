## ui.R

fluidPage(
  pageWithSidebar(
    headerPanel("Analaysis and visualization of JUMP quantification results"),
    sidebarPanel(
      h2("About"),
      p("This tool draws a volcano plot- typically used for showing relationships between fold change and p-values of genes (probes) in microarray studies.",
        "Points are sized according to their X-value, and colour intensity is by their Y-value"), 
      p("Which genes to labeSl are chosen based on Euclidean distance from the origin, but can be biased towards the X or Y dimension. All other settings should be self explanatory.",
        "Editing any settings will cause the plot to be redrawn, label positions are fixed above the point currently."),
      fileInput("inputFile", "Choose File",
                accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      numericInput("nGroups", "Number of groups in samples", value = 2, min = 2)
    ),
      
    mainPanel(
      tabsetPanel(
        tabPanel("Select samples for each group", 
                uiOutput("sampleGroups"), 
                actionButton("submitGroups", "Update")),
        tabPanel("Volcano Plot", plotOutput(outputId = 'volcano', width = "800px", height = "800px")),
        tabPanel("Histogram", plotOutput(outputId = 'histogram', width = "100%")),
        tabPanel("Data Preview", tableOutput("contents"))
      )
    )
  )
)