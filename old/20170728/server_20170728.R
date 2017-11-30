library(shiny)
library(readxl)
library(pheatmap)
library(RColorBrewer)

function (input, output) {
  ## Load jump -q parameter file
  
  
  ## Data frame of differentially expressed peptides/proteins
  df_DE = reactive({
    ## Data loading
    inputFile = switch(input$level,
           "Unique peptides" = "id_uni_pep_quan.xlsx",
           "All peptides" = "id_all_pep_quan.xlsx",
           "Unique proteins" = "id_uni_prot_quan.xlsx",
           "All proteins" = "id_all_prot_quan.xlsx")
    if (grep("pep", inputFile) == 1) {
      df = read_excel(inputFile, skip = 4) # Peptide publication table
    } else {
      df = read_excel(inputFile, skip = 1) # Protein publication table
    }
    colIndsFold = grep('Log2Fold', colnames(df))
    colNamesFold = colnames(df)[colIndsFold]
    isMultipleGroups = length(colIndsFold) > 1

    ## FDR and log2FC threshold
    thrFDR = as.numeric(input$fdr)
    thrFold = as.numeric(input$log2fc)

    ## Selection of differentially expressed peptides/proteins and generate a heatmap
    df = data.frame(df) # change the class from tbl_df to data.frame
    if (isMultipleGroups) {
      logFC = apply(df[, colIndsFold], 1, min) # the minimum fold-change of multiple groups
      rowIndDE = df$FDR < thrFDR & abs(logFC) > thrFold
    } else {
      logFC = as.numeric(df[, colIndsFold])
      rowIndDE = df$FDR < thrFDR & abs(logFC) > thrFold
    }
    colIndDE = grep('sig', colnames(df)) # 1st column = entry (peptide or protein), 2nd column ~ = expression of reporter ions
    df_DE = df[rowIndDE, colIndDE]
    rownames(df_DE) = df[rowIndDE, 1]
    colNamesVec = colnames(df_DE)
    colNamesVec = gsub('\\.\\.', '_', colNamesVec)
    colNamesVec = gsub('\\.', '', colNamesVec)
    colnames(df_DE) = colNamesVec
    df_DE_matrix = t(scale(t(df_DE)))

    bk = c(min(df_DE_matrix), seq(-2, 2, length.out = 99), max(df_DE_matrix))
    col = c("blue", colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(98), "red")
    pheatmap(df_DE_matrix, main = paste("Differentially expressed peptides/proteins,", dim(df_DE_matrix)[1], "entries"),
             treeheight_col = 20, show_rownames = (sum(rowIndDE) < 50),
             breaks = bk, colors = col,
             cellwidth = 30, fontsize = 12)
  })
  output$heatmap = renderPlot(df_DE())
}













