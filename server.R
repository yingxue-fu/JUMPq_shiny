## server.R
library(readxl)

function (input, output) {
  ## Load JUMP -q output file (either id_uni_pep_quan.xlsx or id_uni_prot_quan.xlsx)
  inputData = reactive({
    inFile = input$inputFile
    if (length(grep("pep", inFile))) {
      tbl = read_excel(inFile, skip = 4) # Peptide publication table
    } else {
      tbl = read_excel(inFile, skip = 1) # Protein publication table
    }
    df = as.data.frame(tbl)
  })
  
  ## Specificiation of sample groups
  inputNumberGroups = reactive(as.integer(input$nGroups))
  output$sampleGroups = renderUI({
    nGroups = inputNumberGroups()
    data = inputData()
    colSampleNames = grep('sig', colnames(data))
    sampleNames = colnames(data)[colSampleNames]
    lapply (1:nGroups, function(i) {
      checkboxGroupInput(inputId = paste0("Group", i), label = paste("Group", i),
                         choiceNames = as.list(sampleNames), choiceValues = as.list(sampleNames))
    })
  })
  
  # Comparison between groups specified by the user input
  DE = eventReactive(input$submitGroups, {
    ## Load the data
    data = inputData()
    nGroups = inputNumberGroups()
    
    ## Get the information of samples to be compared
    comparison = as.character()
    for (g in 1:nGroups) {
      groupName = paste0("Group", g)
      comparison[g] = paste(input[[groupName]], collapse = ",")
    }
    groups = list()
    nSamples = 0
    compSamples = NULL
    for (g in 1:nGroups) {
      groups[[g]] = unlist(strsplit(comparison[g], ","))
      nSamples = nSamples + length(groups[[g]])
      compSamples = c(compSamples, groups[[g]])
    }
    
    ## Generate a design matrix (which contains the information of comparison)
    subColNames = colnames(data)[which(colnames(data) %in% compSamples)]
    design = matrix(0, nrow = nSamples, ncol = nGroups)
    for (g in 1:nGroups) {
      design[which(subColNames %in% groups[[g]]), g] = 1
    }
    colnames(design) = paste("group", seq(1, nGroups), sep = "")

    ## Generate a contrast matrix and new column names for the LIMMA result table
    contVec = NULL
    newColNames = NULL
    combMatrix = combn(seq(1, nGroups), 2)
    for (j in 1:ncol(combMatrix)) {
      contVec = c(contVec, paste(paste("group", combMatrix[1, j], sep = ""), paste("group", combMatrix[2, j], sep = ""), sep = "-"))
      newColNames = c(newColNames, paste(comparison[combMatrix[1, j]], "/", comparison[combMatrix[2, j]], sep = ""))
    }
    contMatrix = makeContrasts(contrasts = contVec, levels = design)

  })
  
  # output$selectedSamples = renderText({
  #   groupName = paste0("Group", 1)
  #   paste(input[[groupName]])
  # })
  
  # output$txt <- renderText({
  #   icons <- paste(input$icons, collapse = ", ")
  #   paste("You chose", icons)
  # })
  
  output$selectedSamples = renderPrint({
    # data = DE()
  })
  
}