rm(list = ls())

library(readxl)
inputFile = "id_uni_prot_quan.xlsx"
if (length(grep("pep", inputFile))) {
  tbl = read_excel(inputFile, skip = 4) # Peptide publication table
} else {
  tbl = read_excel(inputFile, skip = 1) # Protein publication table
}






inputComparison = "Synaptosome_vs_Tissue"

colIndSamples = grep("sig", colnames(tbl))
colIndPval = grep(paste("p-value", inputComparison, sep = "_"), colnames(tbl))
colIndFDR = grep(paste("FDR", inputComparison, sep = "_"), colnames(tbl))
colIndsFold = grep('Log2Fold', colnames(tbl))
colIndsFold = colIndsFold[colIndsFold < colIndPval]
isMultipleGroups = length(colIndsFold) > 1

## FDR and log2FC threshold
thrFDR = 0.1
thrFold = 0.7

## Selection of DE peptides/proteins and create a heatmap
df = as.data.frame(tbl)
if (isMultipleGroups) {
  logFC = apply(df[, colIndsFold], 1, min) # the minimum fold-change over multiple comparisons
  pval = as.numeric(df[, colIndPval])
  fdr = as.numeric(df[, colIndFDR])
  rowIndDE = fdr < thrFDR & abs(logFC) > thrFold
} else {
  logFC = as.numeric(df[, colIndsFold])
  pval = as.numeric(df[, colIndPval])
  fdr = as.numeric(df[, colIndFDR])
  rowIndDE = fdr < thrFDR & abs(logFC) > thrFold
}
colIndDE = grep("sig", colnames(df))
df_DE = df[rowIndDE, colIndDE]