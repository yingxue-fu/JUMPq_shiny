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
            p("This panel provides explorative data analysis results for your proteomics dataset."),
            p("1. PCA (principal component analysis) plot of samples"),
            p("2. Hierarchical clustering result of samples and peptides/proteins"),
            br(),
            fileInput("inputFile", label = HTML("Choose a file<br />
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
    
    ## Differential expression analysis tab
    tabPanel(
      "Differential expression analysis",
      br(),
      mainPanel(
        width = 12,
        fluidRow(
          column (3, wellPanel(
            p("This panel provides differential expression analysis results for your proteomics dataset."),
            p("1. For two group comparison, a volcano plot and heatmap will be provided"),
            p("2. For three or more group comparison, a heatmap will be provided"),
            br(),
            fileInput("diffFile", label = HTML("Choose a file<br />
                                              id_uni_pep_quan.xlsx or<br />
                                                id_uni_prot_quan.xlsx"),
                      accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
            numericInput("diffNumGroups", label = "Number of groups", value = 2),
            uiOutput("diffGroups"),
            selectInput("diffMetric", label = "Select the measure of significance", 
                        choice = list("p-value" = "p-value",
                                      "FDR" = "FDR"),
                        selected = 1),
            numericInput("diffSigCutoff", label = "Significance level", min = 0, max = 1, step = 0.01, value = 0.05),
            numericInput("diffLogFC", label = "Log2-fold cutoff", value = 1),
            actionButton("diffSubmit", "Submit")
            )),
          column (5, align = "center", wellPanel(
            # textOutput("pcaText"),
            # tags$head(
            #   tags$style(
            #     "#pcaText{font-size: 20px;}"
            #   )
            # ),
            # br(),
            # plotOutput("pcaPlot", height = "600px")
          )),
          column (4, align = "center", wellPanel(
            # textOutput("hclustText"),
            # tags$head(
            #   tags$style(
            #     "#hclustText{font-size: 20px;}"
            #   )
            # ),
            # br(),
            plotOutput("hclustDE", height = "600px")
          ))
        )
      )
    )
    # tabPanel("Volcano Plot", plotOutput(outputId = 'volcano', width = "800px", height = "800px")),
    # tabPanel("Histogram", plotOutput(outputId = 'histogram', width = "100%")),
    # tabPanel("Data Preview", tableOutput("contents"))
  )
)
