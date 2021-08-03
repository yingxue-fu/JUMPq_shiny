library(shiny)
library(readxl)
library(scatterD3)
library(gplots)
library(ggplot2)
library(plotly)
library(heatmaply)
library(InteractiveComplexHeatmap)
library(ComplexHeatmap)
library(circlize)
library(DT)
library(curl)

source("preprocess.R")
source("statTest.R")
server = function (input, output, session) {
    # Increase the maximum size of uploaded file (up to 500MB)
    options(shiny.maxRequestSize = 500 * (1024 ^ 2))
    
    ####################################################
    # Unsupervised analysis, i.e. explorative analysis #
    ####################################################
    # Load JUMP -q output file (either id_uni_pep_quan.xlsx or id_uni_prot_quan.xlsx)
    data1 = reactive ({
        inFileName = input$inputFile1$name
        if (length(grep("pep", inFileName))) {
            tbl = read_excel(input$inputFile1$datapath, skip = 4) # Peptide publication table
            level = "peptide"
        } else {
            tbl = read_excel(input$inputFile1$datapath, skip = 1) # Protein publication table
            level = "protein"
        }
        list(data = as.data.frame(tbl), level = level)
    })
    
    # Load a file containing sample information
    # Column 1: ID = should be the same as the sample labels in JUMP -q output file
    # Column 2~: Any column name is possible, and sample information can be put as text
    metaData1 = reactive ({
        if (!is.null(input$metaFile1)) {
            read.table(input$metaFile1$datapath, sep="\t", header=T)
        } else {
            NULL
        }
    })
    
    # Selection of a data subset (highly variable), This subset is for analyses
    subData1 = eventReactive(input$submit1, {
        df = data1()$data
        dfSample = metaData1()
        level = data1()$level
        metric = input$metric1
        pct = input$variant1
        res = preprocess(df, level, metric, pct)
        return (list(rawData = res$rawData, data = res$data, sampleInfo = dfSample))
    })
    
    ############
    # PCA plot #
    ############
    output$pcaPlot = renderScatterD3({
        df = subData1()$data
        dfSample = subData1()$sampleInfo
        res = prcomp(t(df), center = TRUE, scale = TRUE)
        eigs = res$sdev ^ 2
        res = data.frame(res$x[, 1:2])
        xlab = paste0("PC1 (", round((eigs[1] / sum(eigs)) * 100, 2),"%)")
        ylab = paste0("PC2 (", round((eigs[2] / sum(eigs)) * 100, 2),"%)")
        
        # Plot options
        if (!is.null(input$pcaPointColor)) {
            col_var = dfSample[[input$pcaPointColor]]
            col_lab = input$pcaPointColor
        } else {
            col_var = NULL
            col_lab = "None"
        }
        point_size = input$pcaPointSize
        labels_size = input$pcaLabelSize
        if (!is.null(input$pcaNoLabel)) {
            if (input$pcaNoLabel == TRUE) {
                labels_size = 0
            }
        }
        point_opacity = input$pcaPointOpacity
        
        scatterD3(x = res$PC1,
                  y = res$PC2,
                  lab = rownames(res),
                  xlab = xlab,
                  ylab = ylab,
                  key_var = rownames(res),
                  col_var = col_var,
                  col_lab = col_lab,
                  zoom_on = NULL,
                  zoom_on_level = 3,
                  labels_positions = "auto",
                  point_opacity = point_opacity,
                  point_size = point_size,
                  labels_size = labels_size,
                  transitions = TRUE,
                  left_margin = 50,
                  lasso = TRUE)
    })
    
    # Control panel of the PCA plot
    observeEvent(input$submit1, {
        output$pcaPointColor = renderUI({
            df = metaData1()
            vars = setdiff(colnames(df), "ID")
            vars = c("None", vars)
            selectInput("pcaPointColor", "Color mapping variable",
                        choices = vars, selected = "None")
        })
        output$pcaPointSize = renderUI({
            sliderInput("pcaPointSize", "Symbol size",
                        min = 10, max = 100, step = 10, value = 50)
        })
        output$pcaLabelSize = renderUI({
            sliderInput("pcaLabelSize", "Label size",
                        min = 10, max = 20, step = 1, value = 12)
        })
        output$pcaHideLabel = renderUI({
            checkboxInput("pcaNoLabel", "Hide labels", FALSE)
        })
        output$pcaOpacity = renderUI({
            sliderInput("pcaPointOpacity", "Point opacity",
                        min = 0, max = 1, step = 0.1, value = 1)
        })
    })

    ##########################
    # Heatmap and dendrogram #
    ##########################
    output$hclustPlot = renderPlot({
        df = subData1()$data
        dfSample = subData1()$sampleInfo
        mat = as.matrix(df)
        mat = t(scale(t(mat), center = T, scale = T))
        
        # Heatmap options
        limVal = round(min(abs(min(mat)), abs(max(mat))))
        myColor = colorRamp2(c(-limVal, 0, limVal), c("blue", "white", "red"), space = "RGB")
        clusterColumns = T
        if (!is.null(input$hClusterColumn)) {
            if (input$hClusterColumn == "No") {
                clusterColumns = F
            }
        }
        topAnnot = NULL
        if (!is.null(input$hClusterTA)) {
            if (input$hClusterTA != "None") {
                topAnnot = HeatmapAnnotation(sample_information = dfSample[[input$hClusterTA]], name = input$hClusterTA)
            }
        }
        clusterDistance = "euclidean"
        if (!is.null(input$hClusterDistance)) {
            clusterDistance = input$hClusterDistance
        }
        clusterMethod = "ward.D2"
        if (!is.null(input$hClusterMethod)) {
            clusterMethod = input$hClusterMethod
        }
        rowSplit = NULL
        if (!is.null(input$hClusterRows)) {
            if (input$hClusterRows > 1) {
                rowSplit = input$hClusterRows
            }
        }
        colSplit = NULL
        if (!is.null(input$hClusterColumns)) {
            if (input$hClusterColumns > 1) {
                colSplit = input$hClusterColumns
            }
        }
        
        ht = Heatmap(mat,
                     col = myColor,
                     cluster_columns = clusterColumns,
                     show_row_names = F,
                     clustering_distance_rows = clusterDistance,
                     clustering_distance_columns = clusterDistance,
                     clustering_method_rows = clusterMethod,
                     clustering_method_columns = clusterMethod,
                     row_dend_width = unit(2, "cm"),
                     border = T,
                     row_split = rowSplit,
                     column_split = colSplit,
                     heatmap_legend_param = list(
                         title = "Scaled intensity",
                         legend_height = unit(5, "cm"),
                         title_position = "leftcenter-rot"
                     ),
                     top_annotation = topAnnot
        )
        ht = draw(ht)
    })
    
    # Control panel of the heatmap
    observeEvent(input$submit1, {
        output$hClusterColumn = renderUI({
            selectInput("hClusterColumn", "Cluster colunms?",
                        choices = c("Yes", "No"), 
                        selected = "Yes")
        })
        output$hClusterTA = renderUI({
            df = metaData1()
            vars = setdiff(colnames(df), "ID")
            vars = c("None", vars)
            selectInput("hClusterTA", "Column annotation variable",
                        choices = vars,
                        selected = "None")
        })
        output$hClusterDistance = renderUI({
            selectInput("hClusterDistance", "Distance metric for clustering",
                        choices = c("euclidean", "manhattan", "pearson", "spearman"),
                        selected = "euclidean")
        })
        output$hClusterMethod = renderUI({
            selectInput("hClusterMethod", "Linkage method for clustering",
                        choices = c("single", "complete", "average", "median", "centroid", "ward.D", "ward.D2"),
                        selected = "ward.D2")
        })
        output$hClusterRows = renderUI({
            sliderInput("hClusterRows", "Split rows (based on dendrogram)",
                        min = 1, max = 10, value = 1)
        })
        output$hClusterColumns = renderUI({
            sliderInput("hClusterColumns", "Split columns (based on dendrogram)",
                        min = 1, max = 10, value = 1)
        })
    })
    
    ##############
    # Data table #
    ##############
    output$dataTable1 = DT::renderDataTable({
        data = subData1()$data
        ## Since data is log2-transformed,
        ## it needs to be re-transformed to raw-scale intensity levels
        ## for showing a data table
        data = round(2 ** data, digits = 2)
    }, selection = 'single', options = list(autoWidth = FALSE, scrollX = TRUE, pageLength = 5))
    
    # Plot of the selected rows from the data table
    output$plotDataTable1 = renderPlot({
        data = subData1()$data
        
        # Since data is log2-transformed, it needs to be re-transformed to the raw-scale intensity levels for a data table
        data = round(2 ** data, digits = 2)
        rowInd = input$dataTable1_rows_selected
        if (length(rowInd) == 1) {
            x = as.numeric(data[rowInd, ])
            df = data.frame(samples = colnames(data), intensity = x)
            g = ggplot(df, aes(x = samples, y = intensity)) +
                geom_bar(stat = "identity") +
                theme(text = element_text(size = 15),
                      axis.text.x = element_text(angle = 90, hjust = 1)) +
                scale_x_discrete(limits = colnames(data)) +
                coord_cartesian(ylim = c(0.8 * min(x), max(x)))
            plot(g)
        } else {
            rowInd = 1
            x = as.numeric(data[rowInd, ])
            df = data.frame(samples = colnames(data), intensity = x)
            g = ggplot(df, aes(x = samples, y = intensity)) +
                geom_bar(stat = "identity") +
                theme(text = element_text(size = 15),
                      axis.text.x = element_text(angle = 90, hjust = 1)) +
                scale_x_discrete(limits = colnames(data)) +
                coord_cartesian(ylim = c(0.8 * min(x), max(x)))
            plot(g)
        }
    })
    
    # Download the subset of data (exploratory analysis)
    output$download1 = downloadHandler(
        filename = "exploratory_subset.txt",
        content = function(file) {
            write.table(subData1()$rawData, file, sep = "\t", row.names = FALSE)
        }
    )
    
    ##############################################################
    # Supervised analysis, i.e. differential expression analysis #
    ##############################################################
    # Load JUMP -q output file (either id_uni_pep_quan.xlsx or id_uni_prot_quan.xlsx)
   
    data2 = reactive ({
        inFileName = input$inputFile2$name
        if (length(grep("pep", inFileName))) {
            tbl = read_excel(input$inputFile2$datapath, skip = 4) # Peptide publication table
            level = "peptide"
        } else {
            tbl = read_excel(input$inputFile2$datapath, skip = 1) # Protein publication table
            level = "protein"
        }
        list(data = as.data.frame(tbl), level = level)
    })
    
    # Specificiation of groups of samples
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
    
    # Differentially expressed peptides/proteins
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
    
    # A subset of data selected based on "statRes"
    subData2 = eventReactive(input$submit2, {
        nGroups = nGroups()
        statRes = statRes()
        
        # Load data (for download and analysis, separately)
        rawData = data2()$data
        data = statRes$data
        
        # Handle threshold inputs
        logFC = input$logfc2
        sigMetric = input$metric2
        sigCutoff = input$cutoff2
        resLogFC = statRes$res[, grep("Log2Fold", colnames(statRes$res))]
        if (nGroups > 2) {
            absLogFC = apply(cbind(abs(apply(resLogFC, 1, min)), abs(apply(resLogFC, 1, max))), 1, max)
        } else {
            absLogFC = abs(resLogFC)
        }
        
        # Select DE peptides/proteins and organize a dataset for subsequent analyses
        rowInd = which(statRes$res[[sigMetric]] < sigCutoff & absLogFC >= logFC)
        if (nGroups == 2) {
            data = cbind(data, `p-value` = statRes$res$`p-value`, FDR = statRes$res$FDR, Log2Fold = resLogFC)
        } else if (nGroups > 2) {
            data = cbind(data, `p-value` = statRes$res$`p-value`, FDR = statRes$res$FDR)
            data = cbind(data, resLogFC)
        }
        data = data[rowInd, ]
        data = data[order(data$`p-value`), ]
        
        # Organize a dataset for downloading
        colInd = max(grep("sig[0-9]{3}", colnames(rawData))) # Last column index for reporters
        rawData = rawData[, 1:colInd] # Remove statistical analysis results from jump -q
        if (nGroups == 2) {
            rawData = cbind(rawData, `p-value` = statRes$res$`p-value`, FDR = statRes$res$FDR, Log2Fold = resLogFC)
        } else if (nGroups > 2) {
            rawData = cbind(rawData, `p-value` = statRes$res$`p-value`, FDR = statRes$res$FDR)
            rawData = cbind(rawData, resLogFC)
        }
        rawData = rawData[rowInd, ]
        return (list(rawData = rawData, data = data))
    })
    
    # Volcano plot of differential expression analysis ("statRes" is directly used)
    output$volcanoPlot = renderPlotly({
        # Preparation of the statistical testing result for visualization
        res = statRes()$res
        logFC = input$logfc2
        sigMetric = input$metric2
        sigCutoff = input$cutoff2
        if (sigMetric == "p-value") {
            # res = res[, 2:3]
            res[, 4] = NULL
            ylab = "-log10(p-value)"
        } else if (sigMetric == "FDR") {
            # res = res[, c(2, 4)]
            res[, 3] = NULL
            ylab = "-log10(FDR)"
        }
        colnames(res) = c("id", "logfc", "significance")
        # res[, 2] = -log10(res[, 2])
        xlab = "log2(fold-change)"
        xmin = min(res$logfc)
        xmax = max(res$logfc)
        ymin = 0
        ymax = max(-log10(res$significance))
        
        # # Parameter setup for ggplot
        # ratioDisplay = 4/3
        # ratioValue = (max(res[, 1]) - min(res[, 1])) / (max(res[, 2]) - min(res[, 2]))
        # v = ggplot(data = res, aes(logfc, significance)) +
        #     geom_point(alpha = 0.2, size = 2) +
        #     geom_hline(aes(yintercept = -log10(sigCutoff))) +
        #     geom_vline(aes(xintercept = -logFC)) +
        #     geom_vline(aes(xintercept = logFC)) +
        #     labs(x = xlab, y = ylab) +
        #     coord_fixed(ratioValue / ratioDisplay) +
        #     theme(text = element_text(size = 20))
        # plot(v)
        
        fig = plot_ly(data = res)
        fig = fig %>% add_trace(x = ~logfc, y = ~significance,
                                text = res[, 1],
                                type = "scatter", 
                                mode = "markers",
                                marker = list(size = 5, color = 'rgb(50, 50, 50)', opacity = 0.5),
                                hovertemplate = paste("ID: %{text}",
                                                      "<br>Log2fold: %{x:.4f}", 
                                                      "<br>sig. level: %{y}<extra></extra>"),
                                showlegend = FALSE
                                )
        fig = fig %>% add_segments(x = -logFC, xend = -logFC, y = ymin, yend = ymax, 
                                   line = list(dash = "dash", color = "black"),
                                   showlegend = FALSE)
        fig = fig %>% add_segments(x = logFC, xend = logFC, y = ymin, yend = ymax, 
                                   line = list(dash = "dash", color = "black"),
                                   showlegend = FALSE)
        fig = fig %>% add_segments(x = xmin, xend = xmax, y = sigCutoff, yend = sigCutoff, 
                                   line = list(dash = "dash", color = "black"),
                                   showlegend = FALSE)
        fig = fig %>% layout(yaxis = list(type = "log", autorange = "reversed", tickformat = ".2e"),
                             showlegend = FALSE)
    })
    
    # Heatmap of differentially expressed peptides/proteins (differentially expressed elements selected by "statRes" should be used)
    output$hclustDE = renderPlot({
        data = subData2()$data
        colInd = grep('^sig[0-9]{3}', colnames(data))
        data = data[, colInd]
        mat = as.matrix(data)
        mat = t(scale(t(mat), center = T, scale = F)) # Only mean-centering
        
        # Heatmap options
        limVal = round(min(abs(min(mat)), abs(max(mat))))
        myColor = colorRamp2(c(-limVal, 0, limVal), c("blue", "white", "red"), space = "RGB")
        clusterColumns = T
        if (!is.null(input$hClusterColumn2)) {
            if (input$hClusterColumn2 == "No") {
                clusterColumns = F
            }
        }
        clusterDistance = "euclidean"
        if (!is.null(input$hClusterDistance2)) {
            clusterDistance = input$hClusterDistance2
        }
        clusterMethod = "ward.D2"
        if (!is.null(input$hClusterMethod2)) {
            clusterMethod = input$hClusterMethod2
        }
        rowSplit = NULL
        if (!is.null(input$hClusterRows2)) {
            if (input$hClusterRows2 > 1) {
                rowSplit = input$hClusterRows2
            }
        }
        colSplit = NULL
        if (!is.null(input$hClusterColumns2)) {
            if (input$hClusterColumns2 > 1) {
                colSplit = input$hClusterColumns2
            }
        }
        
        ht = Heatmap(mat,
                     col = myColor,
                     cluster_columns = clusterColumns,
                     show_row_names = F,
                     clustering_distance_rows = clusterDistance,
                     clustering_distance_columns = clusterDistance,
                     clustering_method_rows = clusterMethod,
                     clustering_method_columns = clusterMethod,
                     row_dend_width = unit(2, "cm"),
                     border = T,
                     row_split = rowSplit,
                     column_split = colSplit,
                     heatmap_legend_param = list(
                         title = "Scaled intensity",
                         legend_height = unit(5, "cm"),
                         title_position = "leftcenter-rot"
                     ),
        )
        ht = draw(ht)
    })
    
    # Control panel of the heatmap
    observeEvent(input$submit2, {
        output$hClusterColumn2 = renderUI({
            selectInput("hClusterColumn2", "Cluster colunms?",
                        choices = c("Yes", "No"), 
                        selected = "Yes")
        })
        output$hClusterDistance2 = renderUI({
            selectInput("hClusterDistance2", "Distance metric for clustering",
                        choices = c("euclidean", "manhattan", "pearson", "spearman"),
                        selected = "euclidean")
        })
        output$hClusterMethod2 = renderUI({
            selectInput("hClusterMethod2", "Linkage method for clustering",
                        choices = c("single", "complete", "average", "median", "centroid", "ward.D", "ward.D2"),
                        selected = "ward.D2")
        })
        output$hClusterRows2 = renderUI({
            sliderInput("hClusterRows2", "Split rows (based on dendrogram)",
                        min = 1, max = 10, value = 1)
        })
        output$hClusterColumns2 = renderUI({
            sliderInput("hClusterColumns2", "Split columns (based on dendrogram)",
                        min = 1, max = 10, value = 1)
        })
        
    })
    
    # Data table
    output$dataTable2 = DT::renderDT({
        ## Since data is log2-transformed, it needs to be re-transformed to raw-scale intensity levels for a data table
        data = subData2()$data
        colInd = grep('^sig[0-9]{3}', colnames(data))
        data[, colInd] = round(2 ** data[, colInd], digits = 2)
        data$`p-value` = format(data$`p-value`, digits = 3)
        data$FDR = format(data$FDR, digits = 3)
        colInd = grep('Log2Fold', colnames(data))
        data[, colInd] = round(data[, colInd], digits = 2)
        data
    }, selection = 'single', options = list(autoWidth = FALSE, scrollX = TRUE, pageLength = 5))
    
    # Plot of the selected rows from the data table
    output$plotDataTable2 = renderPlot({
        data = subData2()$data
        colInd = grep('^sig[0-9]{3}', colnames(data))
        data = data[, colInd]
        data = round(2 ** data, digits = 2)
        rowInd = input$dataTable2_rows_selected
        if (length(rowInd) == 1) {
            x = as.numeric(data[rowInd, ])
            df = data.frame(samples = colnames(data), intensity = x)
            g = ggplot(df, aes(x = samples, y = intensity)) +
                geom_bar(stat = "identity") +
                theme(text = element_text(size = 15),
                      axis.text.x = element_text(angle = 90, hjust = 1)) +
                scale_x_discrete(limits = colnames(data)) +
                coord_cartesian(ylim = c(0.8 * min(x), max(x)))
            plot(g)
        } else {
            rowInd = 1
            x = as.numeric(data[rowInd, ])
            df = data.frame(samples = colnames(data), intensity = x)
            g = ggplot(df, aes(x = samples, y = intensity)) +
                geom_bar(stat = "identity") +
                theme(text = element_text(size = 15),
                      axis.text.x = element_text(angle = 90, hjust = 1)) +
                scale_x_discrete(limits = colnames(data)) +
                coord_cartesian(ylim = c(0.8 * min(x), max(x)))
            plot(g)
        }
    })
    
    # Download the subset of data (differential expression analysis)
    output$download2 = downloadHandler(
        filename = "DE_proteins.txt",
        content = function(file) {
            write.table(subData2()$rawData, file, sep = "\t", row.names = FALSE)
        }
    )
}
