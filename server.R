## server.R
rm(list = ls())
library(readxl)
library(gplots)
library(ggplot2)
library(DT)
library(msigdbr)
library(clusterProfiler)
source("R/statTest.R")

function (input, output) {
    ## Increase the maximum size of uploaded file (up to 100MB)
    options(shiny.maxRequestSize = 100 * (1024 ^ 2))
    ######################################################
    ######################################################
    ## Unsupervised analysis, i.e. explorative analysis ##
    ######################################################
    ######################################################
    ## Load JUMP -q output file (either id_uni_pep_quan.xlsx or id_uni_prot_quan.xlsx)
    data1 = reactive ({
        # ## Desktop version
        # inFileName = input$inputFile1$name
        # if (length(grep("pep", inFileName))) {
        #     tbl = read_excel(inFileName, skip = 4) # Peptide publication table
        #     level = "peptide"
        # } else {
        #     tbl = read_excel(inFileName, skip = 1) # Protein publication table
        #     level = "protein"
        # }
        
        ## Server version
        inFileName = input$inputFile1$name
        if (length(grep("pep", inFileName))) {
            tbl = read_excel(input$inputFile1$datapath, skip = 4) # Peptide publication table
            level = "peptide"
        } else {
            tbl = read_excel(input$inputFile1$datapath, skip = 1) # Protein publication table
            level = "protein"
        }
        
        ## Selection of a subset of data according to input parameters
        list(data = as.data.frame(tbl), level = level)
    })
    
    ##################################################
    ## Selection of a data subset (highly variable) ##
    ## This subset is for analyses                  ##
    ##################################################
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
        ## Organize a dataset for subsequent analyses
        ##   Column 1: entry (either protein accession or peptide sequence)
        ##   Column 2~ : log2-transformed intensity values for reporter ions
        entry = entry[rowInd]
        data = log(rawData[rowInd, colInd], 2)
        rownames(data) = entry
        ## Organize a dataset for downloading
        colInd = max(grep("sig", colnames(rawData)))
        rawData = rawData[, 1:colInd] ## Remove statistical analysis results from previous jump -q
        rawData = rawData[rowInd, ]
        return (list(rawData = rawData, data = data))
    })
    
    # ####################################################################
    # ## Selection of a data subset (highly variable)                   ##
    # ## This subset is for showing a data table and downloading a file ##
    # ####################################################################
    # subData1_for_download = eventReactive(input$submit1, {
    #     rawData = data1()$data
    #     
    #     ## CV or MAD calcuation is based on log2-transformed intensities
    #     ## but output format is raw-intensity scale
    #     colInd = grep("sig", colnames(rawData))
    #     data = log(rawData[, colInd], 2)
    #     cv = apply(data, 1, sd) / rowMeans(data)
    #     mad = apply(abs(data - apply(data, 1, median)), 1, median)
    #     threshold = as.numeric(input$variant1)/ 100 ## Threshold percentage
    #     rowInd = NULL
    #     if (as.numeric(input$metric1) == 1) {
    #         rowInd = cv > quantile(cv, prob = 1 - threshold)
    #     } else if (as.numeric(input$metric1) == 2) {
    #         rowInd = mad > quantile(mad, prob = 1 - threshold)
    #     }
    #     rawData[rowInd, ]
    # })
    
    ##############
    ## PCA plot ##
    ##############
    output$pcaPlot = renderPlot({
        reactive(input$submit1, {
            data = subData1()$data
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
                      theme(text = element_text(size = 12),
                            axis.text = element_text(size = 14),
                            axis.title = element_text(size = 14)))
        })
    })
    
    ############################
    ## Heatmap and dendrogram ##
    ############################
    output$hclustPlot = renderPlot({
        reactive(input$submit1, {
            data = subData1()$data
            mat = as.matrix(data)
            mat = t(scale(t(mat), center = T, scale = F)) # Only mean-centering
            ## heatmap.2
            myColor <- colorRampPalette(c("blue", "white", "red"))(n = 100)
            limVal = round(min(abs(min(mat)), abs(max(mat))))
            myBreaks = seq(-limVal, limVal, length.out = 101)
            par(oma = c(10, 3, 1, 3), mar = c(1, 1, 1, 1))
            h = heatmap.2(x = mat, density.info = "n", trace = "n", labRow = F, col = myColor,
                          hclust = function(x) hclust(x, method = "ward.D2"),
                          lhei = c(1, 6.5), lwid = c(2, 10), breaks = myBreaks,
                          key.par = list(mar= c(5, 0, 0, 0)), key.title = NA,
                          key.xlab = "scaled intensity")
        })
    })
    
    ################
    ## Data table ##
    ################
    output$dataTable1 = DT::renderDataTable({
        data = subData1()$data
        ## Since data is log2-transformed,
        ## it needs to be re-transformed to raw-scale intensity levels
        ## for showing a data table
        data = round(2 ** data, digits = 2)
    }, selection = 'single', options = list(scrollX = TRUE))
    
    ###################################################
    ## Plot of the selected rows from the data table ##
    ###################################################
    output$plotDataTable1 = renderPlot({
        data = subData1()$data
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
    
    ########################################################
    ## Download the subset of data (exploratory analysis) ##
    ########################################################
    output$download1 = downloadHandler(
        filename = "exploratory_subset.txt",
        content = function(file) {
            write.table(subData1()$rawData, file, sep = "\t", row.names = FALSE)
        }
    )
    
    ################################################################
    ################################################################
    ## Supervised analysis, i.e. differential expression analysis ##
    ################################################################
    ################################################################
    ## Load JUMP -q output file (either id_uni_pep_quan.xlsx or id_uni_prot_quan.xlsx)
    data2 = reactive ({
        # ## Desktop version
        # inFileName = input$inputFile2$name
        # if (length(grep("pep", inFileName))) {
        #     tbl = read_excel(inFileName, skip = 4) # Peptide publication table
        #     level = "peptide"
        # } else {
        #     tbl = read_excel(inFileName, skip = 1) # Protein publication table
        #     level = "protein"
        # }
        
        ## Server version
        inFileName = input$inputFile2$name
        if (length(grep("pep", inFileName))) {
            tbl = read_excel(input$inputFile2$datapath, skip = 4) # Peptide publication table
            level = "peptide"
        } else {
            tbl = read_excel(input$inputFile2$datapath, skip = 1) # Protein publication table
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
    
    ################################################
    ## Differentially expressed peptides/proteins ##
    ################################################
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
    
    ##################################################
    ## A subset of data selected based on "statRes" ##
    ##################################################
    subData2 = eventReactive(input$submit2, {
        nGroups = nGroups()
        statRes = statRes()
        ## Load data (for download and analysis, separately)
        rawData = data2()$data
        data = statRes$data
        ## Handle threshold inputs
        logFC = input$logfc2
        sigMetric = input$metric2
        sigCutoff = input$cutoff2
        resLogFC = statRes$res[, grep("Log2Fold", colnames(statRes$res))]
        if (nGroups > 2) {
            absLogFC = apply(cbind(abs(apply(resLogFC, 1, min)), abs(apply(resLogFC, 1, max))), 1, max)
        } else {
            absLogFC = abs(resLogFC)
        }
        ## Select DE peptides/proteins
        rowInd = which(statRes$res[[sigMetric]] < sigCutoff & absLogFC >= logFC)
        ## Organize a dataset for subsequent analyses
        # data = cbind(data, `p-value` = statRes$res$`p-value`, FDR = statRes$res$FDR)
        # data = cbind(data, resLogFC)
        if (nGroups == 2) {
            data = cbind(data, `p-value` = statRes$res$`p-value`, FDR = statRes$res$FDR, Log2Fold = resLogFC)
        } else if (nGroups > 2) {
            data = cbind(data, `p-value` = statRes$res$`p-value`, FDR = statRes$res$FDR)
            data = cbind(data, resLogFC)
        }
        data = data[rowInd, ]
        data = data[order(data$`p-value`), ]
        ## Organize a dataset for downloading
        colInd = max(grep("sig[0-9]{3}", colnames(rawData))) ## Last column index for reporters
        rawData = rawData[, 1:colInd] ## Remove statistical analysis results from jump -q
        # rawData = cbind(rawData, `p-value` = statRes$res$`p-value`, FDR = statRes$res$FDR)
        # rawData = cbind(rawData, resLogFC)
        if (nGroups == 2) {
            rawData = cbind(rawData, `p-value` = statRes$res$`p-value`, FDR = statRes$res$FDR, Log2Fold = resLogFC)
        } else if (nGroups > 2) {
            rawData = cbind(rawData, `p-value` = statRes$res$`p-value`, FDR = statRes$res$FDR)
            rawData = cbind(rawData, resLogFC)
        }
        rawData = rawData[rowInd, ]
        return (list(rawData = rawData, data = data))
    })

    ######################################################
    ## Volcano plot of differential expression analysis ##
    ## - "statRes" can be directly used                 ##
    ######################################################
    output$volcanoPlot = renderPlot({
        reactive(input$submit2, {
            ## Preparation of PCA result for visualization
            res = statRes()$res
            logFC = input$logfc2
            sigMetric = input$metric2
            sigCutoff = input$cutoff2
            if (sigMetric == "p-value") {
                res = res[, 2:3]
                ylab = "-log10(p-value)"
            } else if (sigMetric == "FDR") {
                res = res[, c(2, 4)]
                ylab = "-log10(FDR)"
            }
            colnames(res) = c("logfc", "significance")
            res[, 2] = -log10(res[, 2])
            xlab = "log2(fold-change)"
            
            ## Parameter setup for ggplot
            ratioDisplay = 4/3
            ratioValue = (max(res[, 1]) - min(res[, 1])) / (max(res[, 2]) - min(res[, 2]))
            v = ggplot(data = res, aes(logfc, significance)) +
                geom_point(alpha = 0.2, size = 2) + 
                geom_hline(aes(yintercept = -log10(sigCutoff))) + 
                geom_vline(aes(xintercept = -logFC)) + 
                geom_vline(aes(xintercept = logFC)) + 
                labs(x = xlab, y = ylab) +
                coord_fixed(ratioValue / ratioDisplay) + 
                theme(text = element_text(size = 20))
            plot(v)
        })
    })
    
    ############################################################################
    ## Heatmap of differentially expressed peptides/proteins                  ##
    ## Differentially expressed elements selected by "statRes" should be used ##
    ############################################################################
    output$hclustDE = renderPlot({
        reactive(input$submit2, {
            data = subData2()$data
            colInd = grep('sig[0-9]{3}', colnames(data))
            data = data[, colInd]
            mat = as.matrix(data)
            mat = t(scale(t(mat), center = T, scale = F)) # Only mean-centering
            limVal = round(min(abs(min(mat, na.rm = T)), abs(max(mat, na.rm = T))))
            myBreaks = seq(-limVal, limVal, length.out = 101)
            myColor <- colorRampPalette(c("blue", "white", "red"))(n = 100)
            par(oma = c(10, 3, 1, 3), mar = c(1, 1, 1, 1))
            h = heatmap.2(x = mat, density.info = "n", trace = "n", labRow = F, col = myColor,
                          hclust = function(x) hclust(x, method = "ward.D2"),
                          lhei = c(1, 6.5), lwid = c(2, 10), breaks = myBreaks,
                          key.par = list(mar= c(5, 0, 0, 0)), key.title = NA,
                          key.xlab = "scaled intensity")
        })
    })
    
    ############################################################################
    ## Data table                                                             ##
    ## Differentially expressed elements selected by "statRes" should be used ##
    ############################################################################
    output$dataTable2 = DT::renderDT({
        ## Since data is log2-transformed,
        ## it needs to be re-transformed to raw-scale intensity levels
        ## for showing a data table
        data = subData2()$data
        colInd = grep('sig[0-9]{3}', colnames(data))
        data[, colInd] = round(2 ** data[, colInd], digits = 2)
        data$`p-value` = format(data$`p-value`, digits = 3)
        data$FDR = format(data$FDR, digits = 3)
        colInd = grep('Log2Fold', colnames(data))
        data[, colInd] = round(data[, colInd], digits = 2)
        data
    }, selection = 'single', options = list(autoWidth = FALSE, scrollX = TRUE))
    
    ###################################################
    ## Plot of the selected rows from the data table ##
    ###################################################
    output$plotDataTable2 = renderPlot({
        ## Since data is log2-transformed,
        ## it needs to be re-transformed to raw-scale intensity levels
        ## for showing a data table
        data = subData2()$data
        colInd = grep('sig[0-9]{3}', colnames(data))
        data = data[, colInd]
        data = round(2 ** data, digits = 2)
        rowInd = input$dataTable2_rows_selected
        x = as.numeric(data[rowInd, ])
        if (length(rowInd) == 1) {
            par(mar = c(10, 6, 1, 1), mgp = c(5, 1, 0))
            barplot(x, ylab = "intensity", names.arg = colnames(data), las = 2)
        }
    })
    
    ####################################################################
    ## Download the subset of data (differential expression analysis) ##
    ####################################################################
    output$download2 = downloadHandler(
        filename = "DE_proteins.txt",
        content = function(file) {
            write.table(subData2()$rawData, file, sep = "\t", row.names = FALSE)
        }
    )
    
    ######################
    ## Enrichment study ##
    ######################
    enrichmentRes = reactive({
        if (input$enrichmentSubmit == 0) return ()
        # Database species
        if (input$enrichmentSpecies == 1) {
            orgDb ="org.Hs.eg.db"; orgName = "Homo sapiens"
        } else if (input$enrichmentSpecies == 2) {
            orgDb ="org.Mm.eg.db"; orgName = "Mus musculus"
        } else if (input$enrichmentSpecies == 3) {
            orgDb ="org.Rn.eg.db"; orgName = "Rattus norvegicus"
        } else if (input$enrichmentSpecies == 4) {
            orgDb ="org.Sc.sgd.db"; orgName = "Saccharomyces cerevisiae"
        } else if (input$enrichmentSpecies == 5) {
            orgDb ="org.Dm.eg.db"; orgName = "Drosophila melanogaster"
        } else if (input$enrichmentSpecies == 6) {
            orgDb ="org.Ce.eg.db"; orgName = "Caenorhabditis elegans"
        }
        ## Background genes (optional)
        bgGenes = NULL
        bgFileName = input$enrichmentBgGenes$name
        if (!is.null(bgFileName)) { ## Background prots/genes are specified
            bgProts = readLines(bgFileName) ## Desktop version
            # bgProts = readLines(input$enrichmentBgGenes$datapath) ## Server version
            bgGenes = bitr(bgProts, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = orgDb)
        }
        ## Gene set
        if (input$enrichmentGeneset == 1) {
            gsCat = "H"
        } else {
            gsCat = paste0("C", (as.numeric(input$enrichmentGeneset) - 1))
        }
        m_df = msigdbr(species = orgName, category = gsCat)
        m_t2g = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
        ## Convert query protein accessions to gene symbols
        data = subData2()$data
        prots = NULL
        for (i in 1:dim(data)[1]) {
            prots[i] = unlist(strsplit(rownames(data)[i], '\\|'))[2]
        }
        geneSymbols = bitr(prots, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = orgDb)
        geneSymbols = geneSymbols$SYMBOL
        ## Enrichment analysis using "clusterProfiler"
        if (!is.null(bgGenes)) {
            res = enricher(gene = geneSymbols, TERM2GENE = m_t2g, universe = bgGenes,
                           pvalueCutoff = 1, qvalueCutoff = 1)
        } else {
            res = enricher(gene = geneSymbols, TERM2GENE = m_t2g,
                           pvalueCutoff = 1, qvalueCutoff = 1)
        }
        res = res[, c(3, 4, 5, 6, 9, 8)]
        names(res) = c("GeneRatio", "BgRatio", "pvalue", "FDR", "Count", "Genes")
        return (res)
    })

    output$enrichmentTable = DT::renderDataTable({
        if (input$enrichmentSubmit == 0) return ()
        withProgress(message = NULL, detail = "Enrichment analysis", {
            res = enrichmentRes();
            incProgress(1, detail = paste("Done"))})
        # res$pvalue = as.numeric(formatC(res$pvalue, digits = 3))
        # res$FDR = as.numeric(formatC(res$FDR, digits = 3))
        res$pvalue = format(res$pvalue, digits = 2)
        res$FDR = format(res$FDR, digits = 2)
        res
    }, options = list(scrollX = TRUE))
    
    output$enrichmentDownload = downloadHandler(
        filename = "DE_enrichment.txt",
        content = function(file) {
            data = enrichmentRes()
            data = cbind(Geneset = rownames(data), data)
            rownames(data) = NULL
            write.table(data, file, sep = "\t", quote = F, row.names = F)
        }
    )
}