## server.R
library(readxl)
library(ggplot2)
source("R/statTest.R")

function (input, output) {
  ## Load JUMP -q output file (either id_uni_pep_quan.xlsx or id_uni_prot_quan.xlsx)
  inputData = eventReactive(input$inputFile, {
    inFileName = input$inputFile$name
    if (length(grep("pep", inFileName))) {
      tbl = read_excel(inFileName, skip = 4) # Peptide publication table
      level = "peptide"
    } else {
      tbl = read_excel(inFileName, skip = 1) # Protein publication table
      level = "protein"
    }
    list(data = as.data.frame(tbl), level = level)
  })

  ######################################################
  ## Unsupervised analysis, i.e. explorative analysis ##
  ######################################################
  ## PCA plot
  output$pcaPlot = renderPlot({
    reactive(input$inputFile, {
      ## Data reduction
      rawData = inputData()$data
      level = inputData()$level
      if (level == "peptide") {
        entry = rawData[, 1]
      } else if (level == "protein") {
        entry = rawData[, 2]
      }
      ## Definition of a dataset (only includes TMT reporter intensities)
      ## Hereafter, the data is log2-transformed
      colInd = grep("sig", colnames(rawData))
      data = log(rawData[, colInd], 2)
      rownames(data) = entry
      
      ## Preparation of PCA result for visualization
      resPCA = prcomp(t(data), center = TRUE, scale = TRUE)
      eigs = resPCA$sdev ^ 2
      resPCA = data.frame(resPCA$x[, 1:2])
      
      ## Parameter setup for ggplot
      xlabPCA = paste0("PC1 (", round((eigs[1] / sum(eigs)) * 100, 2),"%)")
      ylabPCA = paste0("PC2 (", round((eigs[2] / sum(eigs)) * 100, 2),"%)")
      ratioDisplay = 4/3
      ratioValue = (max(resPCA$PC1) - min(resPCA$PC1)) / (max(resPCA$PC2) - min(resPCA$PC2))
      print(ggplot(data = resPCA[, 1:2], aes(PC1, PC2)) +
              geom_jitter(size = 3) +
              geom_text(aes(label = rownames(resPCA)), vjust = "inward", hjust = "inward") +
              labs(x = xlabPCA, y = ylabPCA) + 
              coord_fixed(ratioValue / ratioDisplay))
    })
  })

  
  
  
  
  
    
  ## Heatmap of data (with hierarchical clustering)
  

  
  # ################################################################
  # ## Supervised analysis, i.e. differential expression analysis ##
  # ################################################################
  # ## Specificiation of sample groups
  # inputNumberGroups = reactive(as.integer(input$nGroups))
  # observeEvent(input$inputFile, {
  #     output$sampleGroups = renderUI({
  #     nGroups = inputNumberGroups()
  #     data = inputData()
  #     colSampleNames = grep('sig', colnames(data))
  #     sampleNames = colnames(data)[colSampleNames]
  #     lapply (1:nGroups, function(i) {
  #       checkboxGroupInput(inputId = paste0("Group", i), label = paste("Group", i),
  #                          choiceNames = as.list(sampleNames), choiceValues = as.list(sampleNames))
  #     })
  #   })
  # })
  # 
  # ## Differentially expressed peptides/proteins
  # DE = eventReactive(input$submitGroups, {
  #   inputFile = input$inputFile$name
  #   level = NULL
  #   if (length(grep("pep", inputFile))) {
  #     level = "peptide"
  #   } else if (length(grep("pro", inputFile))) {
  #     level = "protein"
  #   }
  #   data = inputData()
  #   nGroups = inputNumberGroups()
  #   comparison = as.character()
  #   compSamples = as.character()
  #   for (g in 1:nGroups) {
  #     groupName = paste0("Group", g)
  #     comparison[g] = paste(input[[groupName]], collapse = ",")
  #   }
  #   groups = list()
  #   compSamples = NULL
  #   for (g in 1:nGroups) {
  #     groups[[g]] = unlist(strsplit(comparison[g], ","))
  #     compSamples = c(compSamples, groups[[g]])
  #   }
  #   res = statTest(data, level, comparison)
  # })
  # 
  # ## Volcano plot
  # 
  # ## Heatmap
  # 
  # ## DEpeptide/protein table
  # output$DE = renderTable({
  #   head(DE())
  # })
  
  # DE = eventReactive(input$submitGroups, {
  #   ## Load the data
  #   data = inputData()
  #   nGroups = inputNumberGroups()
  # 
  #   ## Get the information of samples to be compared
  #   comparison = as.character()
  #   for (g in 1:nGroups) {
  #     groupName = paste0("Group", g)
  #     comparison[g] = paste(input[[groupName]], collapse = ",")
  #   }
  #   groups = list()
  #   nSamples = 0
  #   compSamples = NULL
  #   for (g in 1:nGroups) {
  #     groups[[g]] = unlist(strsplit(comparison[g], ","))
  #     nSamples = nSamples + length(groups[[g]])
  #     compSamples = c(compSamples, groups[[g]])
  #   }
  # 
  #   ## Generate a design matrix (which contains the information of comparison)
  #   subColNames = colnames(data)[which(colnames(data) %in% compSamples)]
  #   design = matrix(0, nrow = nSamples, ncol = nGroups)
  #   for (g in 1:nGroups) {
  #     design[which(subColNames %in% groups[[g]]), g] = 1
  #   }
  #   colnames(design) = paste("group", seq(1, nGroups), sep = "")
  # 
  #   ## Generate a contrast matrix and new column names for the LIMMA result table
  #   contVec = NULL
  #   newColNames = NULL
  #   combMatrix = combn(seq(1, nGroups), 2)
  #   for (j in 1:ncol(combMatrix)) {
  #     contVec = c(contVec, paste(paste("group", combMatrix[1, j], sep = ""), paste("group", combMatrix[2, j], sep = ""), sep = "-"))
  #     newColNames = c(newColNames, paste(comparison[combMatrix[1, j]], "/", comparison[combMatrix[2, j]], sep = ""))
  #   }
  #   contMatrix = makeContrasts(contrasts = contVec, levels = design)
  # })
  
  # output$selectedSamples = renderText({
  #   groupName = paste0("Group", 1)
  #   paste(input[[groupName]])
  # })
  
  # output$txt <- renderText({
  #   icons <- paste(input$icons, collapse = ", ")
  #   paste("You chose", icons)
  # })

  
  output$table = renderTable(head(inputData()))
  
  
}