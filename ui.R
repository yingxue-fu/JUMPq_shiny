library(shiny)
library(scatterD3)
library(plotly)

fluidPage(
    # tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")),
    navbarPage(
        "Quantitative analysis of a TMT dataset (jump -q)",
        #############################
        # Exploratory data analysis #
        #############################
        tabPanel(
            "Exploratory data analysis",
            sidebarPanel(
                tags$head(tags$style(HTML('h5 {margin-bottom:0px; margin-top:0px;}'))), width = 3,
                p("1. PCA (principal component analysis) plot of samples", br(),
                  "2. Hierarchical clustering result of samples and peptides/proteins", br(),
                  "3. Data table of highly variant peptides/proteins"), br(),
                fileInput("inputFile1", label = HTML("Choose a file<h5>id_uni_pep_quan.xlsx or id_uni_prot_quan.xlsx</h5>")),
                fileInput("metaFile1", label = HTML("Choose a file<h5>upload a file containing sample information (tab-delimited txt)</h5>")),
                numericInput("variant1", 
                             label = HTML("Proportion of highly variant elements<h5>(if you choose 10, top 10% of highly variant peptides/proteins will be used)</h5>"),
                             value = 10),
                selectInput("metric1", label = "Select the measure of variation",
                            choice = list("Coefficient of variation (CV)" = 1, "Median absolute deviation (MAD)" = 2), selected = 1),
                actionButton("submit1", "Submit")
            ), # end of sidebarPanel

            # Main panel showing the results of the exploratory data analysis
            mainPanel(
                tabsetPanel(
                    tabPanel("Principal component analysis (PCA)", br(),
                             scatterD3Output("pcaPlot", height = "400px", width = "650px"),
                             br(),
                             column(2, offset = 0, uiOutput("pcaPointColor")),
                             column(2, uiOutput("pcaPointSize")),
                             column(2, uiOutput("pcaLabelSize"), uiOutput("pcaHideLabel")),
                             column(2, uiOutput("pcaOpacity")),
                             br(),
                             column(2, downloadButton("downloadPCA", "Download PC1/2 coordinates"))
                    ),
                    tabPanel("Heatmap of the subset of peptides/proteins", br(), align = "center",
                             column(9, plotOutput("hclustPlot", height = "650px", width = "600px")),
                             column(3, uiOutput("hClusterColumn"),
                                    uiOutput("hClusterTA"),
                                    uiOutput("hClusterDistance"),
                                    uiOutput("hClusterMethod"),
                                    uiOutput("hClusterRows"),
                                    uiOutput("hClusterColumns"))
                    ),
                    tabPanel("Data table", br(),
                             DT::dataTableOutput("dataTable1"), br(),
                             downloadButton("download1", "Download"), br(),
                             plotOutput("plotDataTable1")
                    )
                ) # end of tabsetPanel
            ) # end of mainPanel
        ), # end of tabPanel

        ####################################
        # Differential expression analysis #
        ####################################
        tabPanel(
            "Differential expression",
            # Sidebar panel controlling an exploratory data analysis
            sidebarPanel(
                width = 3,
                conditionalPanel(
                    condition = 'input.DE != "Functional enrichment"',
                    p("1. For two group comparison, a volcano plot and heatmap will be provided", br(),
                      "2. For three or more group comparison, a heatmap will be provided", br(),
                      "3. Data table of differentially expressed peptides/proteins", br(),
                      "4. Functional enrichment study of differentially expressed proteins"),
                    br(),
                    fileInput("inputFile2", label = HTML("Choose a file<h5>id_uni_pep_quan.xlsx or id_uni_prot_quan.xlsx</h5>")),
                    fileInput("metaFile2", label = HTML("Choose a file<h5>upload a file containing sample information (tab-delimited txt)</h5>")),
                    uiOutput("groups2"),
                    selectInput("metric2", label = "Select the measure of significance", choice = list("p-value" = "p-value", "FDR" = "FDR"), selected = 1),
                    numericInput("cutoff2", label = "Significance level", min = 0, max = 1, step = 0.01, value = 0.05),
                    numericInput("logfc2", label = "Log2-fold cutoff", value = 1),
                    actionButton("submit2", "Submit")
                )
            ),

            # Main panel showing the results of the differential expression analysis
            mainPanel(
                tabsetPanel(
                    id = "DE",
                    tabPanel("Volcano plot", br(),
                             plotlyOutput("volcanoPlot", height = "500px", width = "700px")
                    ),
                    tabPanel("Heatmap of differentially expressed peptides/proteins", br(), 
                             align = "center", 
                             column(9, plotOutput("hclustDE", height = "650px", width = "600px")),
                             column(3, uiOutput("hClusterColumn2"),
                                    uiOutput("hClusterDistance2"),
                                    uiOutput("hClusterMethod2"),
                                    uiOutput("hClusterRows2"),
                                    uiOutput("hClusterColumns2"))
                    ),
                    tabPanel("Data table", br(), br(), 
                             DT::dataTableOutput("dataTable2"), br(), 
                             downloadButton("download2", "Download"), br(), 
                             plotOutput("plotDataTable2")
                    )
                ) # endof tabsetPanel
            ) # end of mainPanel
        ) # end of tabPanel
    ) # end of navbarPage
) # end of fluidPage
