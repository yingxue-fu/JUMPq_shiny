## ui.R
fluidPage(
    # h1("Analysis and visualization of JUMP quantification results", align = "center"),
    # headerPanel("Analysis and visualization of JUMP quantification results"),
    navbarPage(
        "Quantitative analysis of a TMT dataset (jump -q)",
        ###############################
        ## Exploratory data analysis ##
        ###############################
        tabPanel(
            "Exploratory data analysis",
            ## Sidebar panel controlling an exploratory data analysis
            sidebarPanel(
                width = 3,
                p("1. PCA (principal component analysis) plot of samples", br(),
                  "2. Hierarchical clustering result of samples and peptides/proteins", br(),
                  "3. Data table of highly variant peptides/proteins"),
                fileInput("inputFile1", 
                          label = HTML("Choose a file<h5>id_uni_pep_quan.xlsx or id_uni_prot_quan.xlsx</h5>")),
                numericInput("variant1",
                             label = HTML("Proportion of highly variant elements<h5>(if you choose 10, top 10% of highly variant peptides/proteins will be used)</h5>"),
                             value = 10),
                selectInput("metric1", label = "Select the measure of variation",
                            choice = list("Coefficient of variation (CV)" = 1, "Median absolute deviation (MAD)" = 2), selected = 1),
                actionButton("submit1", "Submit")
            ), ## sidebarPanel
            
            ## Main panel showing the results of the exploratory data analysis
            # mainPanel(
                tabsetPanel(
                    tabPanel("Principal component analysis (PCA)", br(), plotOutput("pcaPlot", height = "500px")),
                    tabPanel("Heatmap of the subset of peptides/proteins", br(), align = "center", plotOutput("hclustPlot", height = "700px", width = "500px")),
                    tabPanel("Data table", br(), DT::dataTableOutput("dataTable1"), br(), downloadButton("download1", "Download"), br(), plotOutput("plotDataTable1"))
                )
            # ) ## mainPanel
        ), ## tabPanel
        
        ######################################
        ## Differential expression analysis ##
        ######################################
        tabPanel(
            "Differential expression",
            ## Sidebar panel controlling an exploratory data analysis
            sidebarPanel(
                width = 3,
                conditionalPanel(
                    condition = 'input.DE == "Functional enrichment"',
                    p("This panel performs functional enrichment analyses using differentially expressed proteins.", br(),
                      "Genesets are from Molecular Signatures Database (MSigDB) curated by Broad institute"),
                    br(),
                    selectInput(
                        "enrichmentSpecies",
                        label = "1. Select species",
                        choice = list("Homo sapiens" = 1,
                                      "Mus musculus" = 2,
                                      "Rattus norvegicus" = 3,
                                      "Saccharomyces cerevisiae" = 4,
                                      "Drosophila melanogaster" = 5,
                                      "Caenorhabditis elegans" = 6),
                        selected = 1),
                    fileInput(
                        "enrichmentBgGenes",
                        label = HTML("2. Background genes<h5>
                                     (If you want to specify background genes,
                                     please upload a .txt file containing
                                     UniProt accession numbers)</h5>")),
                    selectInput(
                        "enrichmentGeneset",
                        label = HTML("3. Select database category<h5>
                                     (The Molecular Signatures Database, MSigDB,
                                     is a collection of 8 annotated gene sets)</h5>"),
                        choice = list("H: hallmark gene sets" = 1,
                                      "C1: positional gene sets" = 2,
                                      "C2: curated gene sets (pathways)" = 3,
                                      "C3: motif gene sets" = 4,
                                      "C4: computational gene sets" = 5,
                                      "C5: GO gene sets" = 6,
                                      "C6: oncogenic gene sets" = 7,
                                      "C7: immunologic gene sets" = 8),
                        selected = 1),
                    actionButton("enrichmentSubmit", "Submit")
                    ),
                conditionalPanel(
                    condition = 'input.DE != "Functional enrichment"',
                    p("1. For two group comparison, a volcano plot and heatmap will be provided", br(),
                      "2. For three or more group comparison, a heatmap will be provided", br(),
                      "3. Data table of differentially expressed peptides/proteins", br(),
                      "4. Functional enrichment study of differentially expressed proteins"),
                    fileInput("inputFile2", label = HTML("Choose a file<h5>id_uni_pep_quan.xlsx or id_uni_prot_quan.xlsx</h5>")),
                    numericInput("numGroups2", label = "Number of groups", value = 2),
                    uiOutput("groups2"),
                    selectInput("metric2", label = "Select the measure of significance", choice = list("p-value" = "p-value", "FDR" = "FDR"), selected = 1),
                    numericInput("cutoff2", label = "Significance level", min = 0, max = 1, step = 0.01, value = 0.05),
                    numericInput("logfc2", label = "Log2-fold cutoff", value = 1),
                    actionButton("submit2", "Submit")
                )
                ),
            
            ## Main panel showing the results of the differential expression analysis
            mainPanel(
                tabsetPanel(
                    id = "DE",
                    tabPanel(
                        "Volcano plot", br(), 
                        conditionalPanel("input.numGroups2 == 2", plotOutput("volcanoPlot")),
                        conditionalPanel("input.numGroups2 > 2", h5("Volcano plot is not available for more than two groups"))
                    ),
                    tabPanel("Heatmap of differentially expressed peptides/proteins", br(), align = "center", plotOutput("hclustDE", height = "700px", width = "500px")),
                    tabPanel("Data table", br(), br(), DT::dataTableOutput("dataTable2"), br(), downloadButton("download2", "Download"), br(), plotOutput("plotDataTable2")),
                    tabPanel("Functional enrichment", br(), br(), DT::dataTableOutput("enrichmentTable"), br(), downloadButton("enrichmentDownload", "Download"))
                )
            )
        )
    )
)
