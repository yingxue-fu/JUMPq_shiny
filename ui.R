## ui.R
fluidPage(
  headerPanel("Analysis and visualization of JUMP quantification results"),
  sidebarPanel(
    h2("About"),
    p("This tool provides exploratory data analysis and differential expression analysis results of a  proteomics dataset."),
    p("Which genes to labeSl are chosen based on Euclidean distance from the origin, but can be biased towards the X or Y dimension. All other settings should be self explanatory.",
      "Editing any settings will cause the plot to be redrawn, label positions are fixed above the point currently."),
    fileInput("inputFile", "Choose File",
              accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
    # numericInput("nGroups", "Number of groups in samples", value = 2, min = 2),
    width = 3
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel("Exploratory analysis",
               sidebarPanel(
                 p("This panel provides explorative data analysis results for your proteomics dataset in three ways."),
                 p("1. PCA (principal component analysis) plot of samples"),
                 p("2. Hierarchical clustering result of samples and peptides/proteins"),
                 p("3. MA plot showing the similarity between any two samples")
               ),
               mainPanel(plotOutput("pcaPlot"))),
      tabPanel("Select samples for each group",
               uiOutput("sampleGroups"),
               actionButton("submitGroups", "Update"),
               uiOutput("DE")),
      tabPanel("Volcano Plot", plotOutput(outputId = 'volcano', width = "800px", height = "800px")),
      tabPanel("Histogram", plotOutput(outputId = 'histogram', width = "100%")),
      tabPanel("Data Preview", tableOutput("contents"))
    )
  )
)