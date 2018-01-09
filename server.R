## server.R
rm(list = ls())
library(readxl)
library(ggplot2)
library(pheatmap)
library(DT)
source("R/statTest.R")

function (input, output) {
  ######################################################
  ## Unsupervised analysis, i.e. explorative analysis ##
  ######################################################
  ## Load JUMP -q output file (either id_uni_pep_quan.xlsx or id_uni_prot_quan.xlsx)
  data1 = reactive ({
    inFileName = input$inputFile1$name
    if (length(grep("pep", inFileName))) {
      tbl = read_excel(inFileName, skip = 4) # Peptide publication table
      level = "peptide"
    } else {
      tbl = read_excel(inFileName, skip = 1) # Protein publication table
      level = "protein"
    }

    ## Selection of a subset of data according to input parameters
    list(data = as.data.frame(tbl), level = level)
  })

  ## Selection of a data subset (highly variable)
  ## This subset is for analyses
  subData1 = eventReactive(input$submit1, {
    rawData = data1()$data
    level = data1()$level
    if (level == "peptide") {
      entry = rawData[, 1]
    } else if (level == "protein") {
      entry = rawData[, 2]
    }

    ## CV or MAD calcuation is based on log2-transformed intensities
    ## but output format is raw-intensity scale
    colInd = grep("sig", colnames(rawData))
    data = log(rawData[, colInd], 2)
    cv = apply(data, 1, sd) / rowMeans(data)
    mad = apply(abs(data - apply(data, 1, median)), 1, median)
    threshold = as.numeric(input$variant1)/ 100 ## Threshold percentage
    rowInd = NULL
    if (as.numeric(input$metric1) == 1) {
      rowInd = cv > quantile(cv, prob = 1 - threshold)
    } else if (as.numeric(input$metric1) == 2) {
      rowInd = mad > quantile(mad, prob = 1 - threshold)
    }

    ## Return data for the following analyses
    ## Column 1: entry (either protein accession or peptide sequence)
    ## Column 2~ : log2-transformed intensity values for reporter ions
    entry = entry[rowInd]
    data = log(rawData[rowInd, colInd], 2)
    rownames(data) = entry
    return (data)
  })


  ## Selection of a data subset (highly variable)
  ## This subset is for showing a data table and downloading a file
  subData1_for_download = eventReactive(input$submit1, {
    rawData = data1()$data

    ## CV or MAD calcuation is based on log2-transformed intensities
    ## but output format is raw-intensity scale
    colInd = grep("sig", colnames(rawData))
    data = log(rawData[, colInd], 2)
    cv = apply(data, 1, sd) / rowMeans(data)
    mad = apply(abs(data - apply(data, 1, median)), 1, median)
    threshold = as.numeric(input$variant1)/ 100 ## Threshold percentage
    rowInd = NULL
    if (as.numeric(input$metric1) == 1) {
      rowInd = cv > quantile(cv, prob = 1 - threshold)
    } else if (as.numeric(input$metric1) == 2) {
      rowInd = mad > quantile(mad, prob = 1 - threshold)
    }
    rawData[rowInd, ]
  })

  ## PCA plot
  output$pcaPlot = renderPlot({
    reactive(input$submit1, {
      data = subData1()
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
              geom_text(aes(label = rownames(resPCA)), vjust = "inward", hjust = "inward", size = 5) +
              labs(x = xlabPCA, y = ylabPCA) +
              coord_fixed(ratioValue / ratioDisplay) +
              theme(text = element_text(size = 12)))
    })
  })

  ## Hierarchical clustering of samples (heatmap)
  output$hclustPlot = renderPlot({
    reactive(input$submit1, {
      data = subData1()
      color <- colorRampPalette(c("blue", "white", "red"))(n = 100)
      mat = as.matrix(data)
      mat = t(scale(t(mat), center = T, scale = F)) # Only mean-centering
      limVal = round(min(abs(min(mat)), abs(max(mat))))
      matBreaks = seq(-limVal, limVal, length.out = 100)
      pheatmap(mat = mat, color = color, breaks = matBreaks, show_rownames = F,
               clustering_method = "ward.D2", fontsize = 12)
    })
  })

  ## Data table
  output$dataTable1 = DT::renderDataTable({
    data = subData1()
    ## Since data is log2-transformed,
    ## it needs to be re-transformed to raw-scale intensity levels
    ## for showing a data table
    data = round(2 ** data, digits = 2)
  }, options = list(scrollX = TRUE))

  ## Plot of the selected rows from the data table
  output$plotDataTable1 = renderPlot({
    data = subData1()
    ## Since data is log2-transformed,
    ## it needs to be re-transformed to raw-scale intensity levels
    ## for showing a data table
    data = round(2 ** data, digits = 2)
    rowInd = input$dataTable1_rows_selected
    x = as.numeric(data[rowInd, ])
    if (length(rowInd) == 1) {
      par(mar = c(10, 6, 1, 1), mgp = c(5, 1, 0))
      barplot(x, ylab = "intensity", names.arg = colnames(data), las = 2)
    }
  })
  
  ## Download the subset of data (exploratory analysis)
  output$download1= downloadHandler(
    filename = "exploratory_subset.txt",
    content = function(file) {
      write.table(subData1_for_download(), file, sep = "\t", row.names = FALSE)
    }
  )

  ################################################################
  ## Supervised analysis, i.e. differential expression analysis ##
  ################################################################
  ## Load JUMP -q output file (either id_uni_pep_quan.xlsx or id_uni_prot_quan.xlsx)
  data2 = reactive ({
    inFileName = input$inputFile2$name
    if (length(grep("pep", inFileName))) {
      tbl = read_excel(inFileName, skip = 4) # Peptide publication table
      level = "peptide"
    } else {
      tbl = read_excel(inFileName, skip = 1) # Protein publication table
      level = "protein"
    }

    ## Selection of a subset of data according to input parameters
    list(data = as.data.frame(tbl), level = level)
  })

  ## Specificiation of groups of samples
  nGroups = reactive(as.integer(input$numGroups2))
  observeEvent(input$inputFile2, {
      output$groups2 = renderUI({
        data = data2()$data
        nGroups = nGroups()
        colSampleNames = grep('sig', colnames(data))
        sampleNames = colnames(data)[colSampleNames]
        lapply (1:nGroups, function(i) {
          checkboxGroupInput(inputId = paste0("Group", i), label = paste("Group", i),
                             choiceNames = as.list(sampleNames), choiceValues = as.list(sampleNames))
      })
    })
  })

  ## Differentially expressed peptides/proteins
  statRes = eventReactive(input$submit2, {
    data = data2()$data
    level = data2()$level
    nGroups = nGroups()
    comparison = as.character()
    compSamples = as.character()
    for (g in 1:nGroups) {
      groupName = paste0("Group", g)
      comparison[g] = paste(input[[groupName]], collapse = ",")
    }
    groups = list()
    compSamples = NULL
    for (g in 1:nGroups) {
      groups[[g]] = unlist(strsplit(comparison[g], ","))
      compSamples = c(compSamples, groups[[g]])
    }
    statTest(data, level, comparison)
  })

  ## A subset of data selected based on "statRes"
  subData2 = eventReactive(input$submit2, {
    nGroups = nGroups()
    statRes = statRes()
    data = statRes$data
    logFC = input$logfc2
    sigMetric = input$metric2
    sigCutoff = input$cutoff2
    resLogFC = statRes$res[, grep("Log2Fold", colnames(statRes$res))]
    if (nGroups > 2) {
      resLogFC = apply(cbind(abs(apply(resLogFC, 1, min)), abs(apply(resLogFC, 1, max))), 1, max)
    } else {
      resLogFC = abs(resLogFC)
    }
    rowInd = which(statRes$res[[sigMetric]] < sigCutoff & resLogFC >= logFC)
    data = data[rowInd, ]
  })
  
  ## Volcano plot of differential expression analysis
  ## - "statRes" can be directly used
  
  ## Heatmap of differentially expressed peptides/proteins
  ## Differentially expressed elements selected by "statRes" should be used
  output$hclustDE = renderPlot({
    reactive(input$submit2, {
      data = subData2()
      mat = as.matrix(data)
      mat = t(scale(t(mat), center = T, scale = F)) # Only mean-centering
      limVal = round(min(abs(min(mat, na.rm = T)), abs(max(mat, na.rm = T))))
      matBreaks = seq(-limVal, limVal, length.out = 100)
      color <- colorRampPalette(c("blue", "white", "red"))(n = 100)
      pheatmap(mat = mat, 
               breaks = matBreaks, show_rownames = F,
               clustering_method = "ward.D2", fontsize = 12)
    })
  })
  
  ## Data table
  ## Differentially expressed elements selected by "statRes" should be used
  output$dataTable2 = DT::renderDataTable({
    ## Since data is log2-transformed,
    ## it needs to be re-transformed to raw-scale intensity levels
    ## for showing a data table
    data = subData2()
    data = round(2 ** data, digits = 2)
  }, options = list(scrollX = TRUE))
  
  ## Plot of the selected rows from the data table
  output$plotDataTable2 = renderPlot({
    ## Since data is log2-transformed,
    ## it needs to be re-transformed to raw-scale intensity levels
    ## for showing a data table
    data = subData2()
    data = round(2 ** data, digits = 2)
    rowInd = input$dataTable2_rows_selected
    x = as.numeric(data[rowInd, ])
    if (length(rowInd) == 1) {
      par(mar = c(10, 6, 1, 1), mgp = c(5, 1, 0))
      barplot(x, ylab = "intensity", names.arg = colnames(data), las = 2)
    }
  })
  
  ## Download the subset of data (exploratory analysis)
  output$download2= downloadHandler(
    filename = "differentially_expressed_subset.txt",
    content = function(file) {
      write.table(subData2_for_table(), file, sep = "\t", row.names = FALSE)
    }
  )
  
  
  # ## Volcano plot
  # 
  # ## Heatmap
  # 
  # ## DEpeptide/protein table
  # output$DE = renderTable({
  #   head(DE())
  # })
  # 
  # DE = eventReactive(input$submitGroups, {
  #   ## Load the data
  #   data = inputData()$data
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

  
  # output$table = renderTable(head(inputData()))
  
  
}