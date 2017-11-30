rm(list = ls())

inputFile = "id_uni_pep_quan.xlsx"
if (length(grep("pep", inputFile))) {
  df = read_excel(inputFile, skip = 4) # Peptide publication table
} else {
  df = read_excel(inputFile, skip = 1) # Protein publication table
}
inputComparison = "twoGroups_S1,S2,S3,S4,S5:S6,S7,S8,S9,S10"
# inputComparison = "threeGroups_S1,S2,S3:S4,S5,S6,S7:S8,S9,S10"

colIndPval = grep(paste("p-value", inputComparison, sep = "_"), colnames(df))
colIndFDR = grep(paste("FDR", inputComparison, sep = "_"), colnames(df))
colIndsFold = grep('Log2Fold', colnames(df))
colIndsFold = colIndsFold[colIndsFold < colIndPval]
isMultipleGroups = length(colIndsFold) > 1

## FDR and log2FC threshold
thrFDR = 0.9
thrFold = 0.6

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





















# ## Data loading
# library(readxl)
# df = read_excel("id_uni_pep_quan.xlsx", skip = 4)
# colIndsFold = grep('Log2Fold', colnames(df))
# colNamesFold = colnames(df)[colIndsFold]
# isMultipleGroups = length(colIndsFold)
# 
# ## FDR and fold-difference threshold
# thrFDR = 0.1
# thrFold = 1
# 
# ## Selection of differentially expressed peptides/proteins and generate a heatmap
# df = data.frame(df) # change the class from tbl_df to data.frame
# if (isMultipleGroups) {
#   logFC = apply(df[, colIndsFold], 1, min) # the minimum fold-change of multiple groups
#   rowIndDE = df$FDR < thrFDR & abs(logFC) > thrFold
# } else {
#   logFC = as.numeric(df[, colIndsFold])
#   rowIndDE = df$FDR < thrFDR & abs(logFC) > thrFold
# }
# colIndDE = grep('sig', colnames(df)) # 1st column = entry (peptide or protein), 2nd column ~ = expression of reporter ions
# df_DE = df[rowIndDE, colIndDE]
# rownames(df_DE) = df[rowIndDE, 1]

# library(pheatmap)
# pheatmap(t(scale(t(df_DE))), treeheight_col = 20, show_rownames = F)

# library(d3heatmap)
# hclustFun = function(x) hclust(x, method = "ward.D2")
# distFun = function(x) dist(x, method = "euclidean")
# d3heatmap(t(scale(t(df_DE))), hclustfun = hclustFun, distfun = distFun,
#           labRow = NULL, colors = "Blues")

# library(gplots)
# hclustFun = function(x) hclust(x, method = "ward.D2")
# distFun = function(x) dist(x, method = "euclidean")
# distFun = function(x) as.dist(1 - cor(t(x)))
# heatmap.2(as.matrix(t(scale(t(df_DE)))), trace = "none", scale = 'none',
#           hclust = hclustFun, dist = distFun,
#           col = redblue(128), labRow = F,
#           key = TRUE, keysize = 1.5)

# library(gplots)
# data(mtcars)
# x = as.matrix(mtcars)
# rc = rainbow(nrow(x), start = 0, end = 0.3)
# cc = rainbow(ncol(x), start = 0, end = 0.3)
