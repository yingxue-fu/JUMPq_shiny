library(shiny)
library(pheatmap)
library(RColorBrewer)
library(readxl)

function (input, output) {
  ## Receive the input level (unique/all peptides/proteins),
  ## and select comparison study
  compNames = reactive({
    inputFile = switch(input$level,
                       "Peptides" = "id_uni_pep_quan.xlsx",
                       "Proteins" = "id_uni_prot_quan.xlsx")
    if (length(grep("pep", inputFile))) {
      df = read_excel(inputFile, skip = 4) # Peptide publication table
    } else {
      df = read_excel(inputFile, skip = 1) # Protein publication table
    }
    compColInds = grep('p-value', colnames(df))
    nComps = length(compColInds)
    compNames = NULL
    for (i in 1:nComps) {
      compNames[i] = sub("p-value_", "", colnames(df)[compColInds[i]])
    }
    as.character(compNames)
  })
  output$compSelection = renderUI({selectInput("comparison",
                                               label = "Choose comparison study",
                                               choices = compNames())})
  
  ## Select differentially expressed peptides/proteins according to input information
  heatmapDE = reactive({
    ## Data loading
    inputFile = switch(input$level,
                       "Peptides" = "id_uni_pep_quan.xlsx",
                       "Proteins" = "id_uni_prot_quan.xlsx")
    if (length(grep("pep", inputFile))) {
      df = read_excel(inputFile, skip = 4) # Peptide publication table
    } else {
      df = read_excel(inputFile, skip = 1) # Protein publication table
    }
    inputComparison = as.character(input$compSelection)
    colIndPval = grep(paste("p-value", inputComparison, sep = "_"), colnames(df))
    colIndFDR = grep(paste("FDR", inputComparison, sep = "_"), colnames(df))
    colIndsFold = grep('Log2Fold', colnames(df))
    colIndsFold = colIndsFold[colIndsFold < colIndPval]
    isMultipleGroups = length(colIndsFold) > 1
    
    ## FDR and log2FC threshold
    thrFDR = as.numeric(input$fdr)
    thrFold = as.numeric(input$log2fc)
    
    ## Selection of DE peptides/proteins and create a heatmap
    df = data.frame(df)
    if (isMultipleGroups) {
      logFC = apply(df[, colIndsFold], 1, min) # the minimum fold-change over multiple comparisons
      pval = as.numeric(df[, colIndPval])
      fdr = as.numeric(df[, colIndFDR])
      rowIndDE = pval < thrFDR & abs(logFC) > thrFold
    } else {
      logFC = as.numeric(df[, colIndsFold])
      pval = as.numeric(df[, colIndPval])
      fdr = as.numeric(df[, colIndFDR])
      rowIndDE = pval < thrFDR & abs(logFC) > thrFold
    }
    colIndDE = grep("sig", colnames(df))
    df_DE = df[rowIndDE, colIndDE]
    df_DE_matrix = t(scale(t(df_DE)))
    # attr(df_DE, "scaled:center") = NULL
    # attr(df_DE, "scaled:scale") = NULL
    pheatmap(df_DE_matrix, treeheight_col = 20, show_rownames = F)
  })
  # output$heatmap = renderDataTable(summary(df_DE()))
  output$heatmap = renderPlot(heatmapDE())
}










