## ui.R
fluidPage(
  h1("Analysis and visualization of JUMP quantification results", align = "center"),
  # headerPanel("Analysis and visualization of JUMP quantification results"),
  tabsetPanel(
    ## Exploratory data analysis tab
    tabPanel(
      "Exploratory data anaysis",
      br(),
      mainPanel(
        width = 12,
        fluidRow(
          column (3, wellPanel(
            p("This panel provides explorative data analysis results for your proteomics dataset in three ways."),
            p("1. PCA (principal component analysis) plot of samples"),
            p("2. Hierarchical clustering result of samples and peptides/proteins"),
            br(),
            fileInput("expFile", label = HTML("Choose a file<br />
                                              id_uni_pep_quan.xlsx or<br />
                                              id_uni_prot_quan.xlsx"),
                      accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
            numericInput("expVariant", label = HTML("Proportion of highly variant elements<br />
                                                    (e.g. if you choose 10, top 10% of highly variant peptides/proteins 
                                                    will be used)"), value = 10),
            selectInput("expMetric", label = "Select the measure of variation", 
                        choice = list("Coefficient of variation (CV)" = 1,
                                      "Median absolute deviation (MAD)" = 2),
                        selected = 1),
            actionButton("expSubmit", "Submit")
            )),
          column (5, align = "center", wellPanel(
            textOutput("pcaText"),
            tags$head(
              tags$style(
                "#pcaText{font-size: 20px;}"
              )
            ),
            br(),
            plotOutput("pcaPlot", height = "600px")
          )),
          column (4, align = "center", wellPanel(
            textOutput("hclustText"),
            tags$head(
              tags$style(
                "#hclustText{font-size: 20px;}"
              )
            ),
            br(),
            plotOutput("hclustPlot", height = "600px")
          ))
        )
      )
    ),
    tabPanel("Select samples for each group",
             uiOutput("sampleGroups"),
             actionButton("submitGroups", "Update"),
             uiOutput("DE")),
    tabPanel("Volcano Plot", plotOutput(outputId = 'volcano', width = "800px", height = "800px")),
    tabPanel("Histogram", plotOutput(outputId = 'histogram', width = "100%")),
    tabPanel("Data Preview", tableOutput("contents"))
  )
)
