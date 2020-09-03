###################################################################
# This version of DESeq2 script is edited by Dr. Bin Zhu at Virginia Commonwealth University
# Email: binzhu0824@gmail.com
###################################################################
###################################################################
##Start##

#load DESeq package;
library("DESeq2")
library("vsn")
library("hexbin")
library("pheatmap")  # using R script 'hclust' to test distance 
library("RColorBrewer")
library("ggplot2")
library("ggrepel")

# set files path
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# Read in your counts table;
cts= read.csv("data.csv", header=TRUE, row.names=1)
design = read.csv("design.csv", header=TRUE, row.names=1)



factor_number <- ncol(design)

if (factor_number == 1) {
  ################ DESeq2 ##################
  # change column names
  col_name <- colnames(design)
  factor_1 = col_name[1]
  
  colnames(design) <- c("factor_1")
  
  design$factor_1 <- as.factor(design$factor_1)
  levels(design$factor_1)  
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = design,
                                design= ~ factor_1)  
  ref1 = toString(design$factor_1[1])
  
  dds$factor_1 <- relevel(dds$factor_1, ref = ref1)   
  
  dds <- DESeq(dds)
  
  
  res=resultsNames(dds) # lists the coefficients
  
  res1 <- gsub("factor_1", factor_1, res)
  res1
  
  n=length(res)
  
  lable_list = as.matrix(read.csv('lable.csv', header = F))
  for (i in 2:n)   {
    res2 <- results(dds, name=res[i])
    write.csv(res2,file = paste(res1[i],".csv"))
    
    # volcano plot
    if (i != 4) {
      volvano_plot <- as.data.frame(res2@listData)
      
      volvano_plot$gene_id = row.names(cts)
      volvano_plot$lable_name = NA
      volvano_plot$lable_name[volvano_plot$gene_id %in% lable_list] <- volvano_plot$gene_id[volvano_plot$gene_id %in% lable_list]
      
      volvano_plot$color <- 'grey'
      volvano_plot$color[volvano_plot$log2FoldChange >= 1 & volvano_plot$padj <= 0.05] = 'green'
      volvano_plot$color[volvano_plot$log2FoldChange <= -1 & volvano_plot$padj <= 0.05] = 'red'
      
      
      volvano_plot$padj <- -log10(volvano_plot$padj)
      
      p=ggplot(data=volvano_plot, aes(x=log2FoldChange, y=padj, color = color, label=lable_name)) + 
        geom_point() + 
        scale_color_manual(values=c("green", "grey", "red")) +
        theme(legend.position = "none") +
        geom_text_repel(aes(label = lable_name),size = 3, colour = 'black') 
      
      print(p)
      
      x_limit_neg <- min(quantile(volvano_plot$log2FoldChange,.01, na.rm =T),-2)
      x_limit_pos <- max(quantile(volvano_plot$log2FoldChange,.99, na.rm =T),2)
      volvano_plot$log2FoldChange[volvano_plot$log2FoldChange < x_limit_neg] = x_limit_neg
      volvano_plot$log2FoldChange[volvano_plot$log2FoldChange > x_limit_pos] = x_limit_pos
      
      y_limit_pos <- max(quantile(volvano_plot$padj,.99, na.rm =T),2)
      volvano_plot$padj[volvano_plot$padj > y_limit_pos] = y_limit_pos
      
      p=ggplot(data=volvano_plot, aes(x=log2FoldChange, y=padj, color = color, label=lable_name)) + 
        geom_point() + 
        scale_color_manual(values=c("green", "grey", "red")) +
        theme(legend.position = "none") +
        geom_text_repel(aes(label = lable_name),size = 3, colour = 'black') 
      
      print(p)
    }
  }
  

  
  ############### Data transformations and visualization ###############
  # transfer data using a variance stabilizing transformation
  vsd <- varianceStabilizingTransformation(dds)  
  
  # this gives log2(n + 1)
  ntd <- normTransform(dds)
  
  ############### Heatmap of the count matrix ###############
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:20]
  df <- as.data.frame(colData(dds)[1])
  colnames(df) <- c(factor_1)
  pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
           cluster_cols=TRUE, annotation_col=df)
  
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:nrow(cts)]
  df <- as.data.frame(colData(dds)[1])
  colnames(df) <- c(factor_1)
  pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
           cluster_cols=TRUE, annotation_col=df)
  
  
  
  ############### Heatmap of the sample-to-sample distances ###############
  sampleDists <- dist(t(assay(vsd)))
  
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$factor_1)     
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  
  
  ############### pincipal component plot of the samples ###############
  pcaData <- plotPCA(vsd, intgroup=c("factor_1"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  #pcaData <- avgdist(cts, sample = )
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  ggplot(pcaData, aes(PC1, PC2, color=factor_1)) +
    geom_point(size=3) +
    labs(color=factor_1) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
  
  
  
} else if (factor_number == 2) {
  # change column names
  col_name <- colnames(design)
  factor_1 = col_name[1]
  factor_2 = col_name[2]
  
  colnames(design) <- c("factor_1","factor_2")
  
  design$factor_1 <- as.factor(design$factor_1)
  design$factor_2 <- as.factor(design$factor_2)
  levels(design$factor_1)                                       
  levels(design$factor_2)                                           
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = design,
                                design= ~ factor_1 + factor_2 + factor_2:factor_1 )  
  ref1 = toString(design$factor_1[1])
  ref2 = toString(design$factor_2[1])
  
  dds$factor_1 <- relevel(dds$factor_1, ref = ref1)   
  dds$factor_2 <- relevel(dds$factor_2, ref = ref2)   
  
  dds <- DESeq(dds)
  
  
  res=resultsNames(dds) # lists the coefficients
  
  res1 <- gsub("factor_1", factor_1, res)
  res1 <- gsub("factor_2", factor_2, res1)
  res1
  
  n=length(res)
  
  lable_list = as.matrix(read.csv('lable.csv', header = F))
  for (i in 2:n)   {
    res2 <- results(dds, name=res[i])
    write.csv(res2,file = paste(res1[i],".csv"))
    
    # volcano plot
    if (i != 4) {
      volvano_plot <- as.data.frame(res2@listData)
      
      volvano_plot$gene_id = row.names(cts)
      volvano_plot$lable_name = NA
      volvano_plot$lable_name[volvano_plot$gene_id %in% lable_list] <- volvano_plot$gene_id[volvano_plot$gene_id %in% lable_list]
      
      volvano_plot$color <- 'grey'
      volvano_plot$color[volvano_plot$log2FoldChange >= 1 & volvano_plot$padj <= 0.05] = 'green'
      volvano_plot$color[volvano_plot$log2FoldChange <= -1 & volvano_plot$padj <= 0.05] = 'red'
      
      
      volvano_plot$padj <- -log10(volvano_plot$padj)
      
      p=ggplot(data=volvano_plot, aes(x=log2FoldChange, y=padj, color = color, label=lable_name)) + 
        geom_point() + 
        scale_color_manual(values=c("green", "grey", "red")) +
        theme(legend.position = "none") +
        geom_text_repel(aes(label = lable_name),size = 3, colour = 'black') 
      
      print(p)
      
      x_limit_neg <- min(quantile(volvano_plot$log2FoldChange,.01, na.rm =T),-2)
      x_limit_pos <- max(quantile(volvano_plot$log2FoldChange,.99, na.rm =T),2)
      volvano_plot$log2FoldChange[volvano_plot$log2FoldChange < x_limit_neg] = x_limit_neg
      volvano_plot$log2FoldChange[volvano_plot$log2FoldChange > x_limit_pos] = x_limit_pos
      
      y_limit_pos <- max(quantile(volvano_plot$padj,.99, na.rm =T),2)
      volvano_plot$padj[volvano_plot$padj > y_limit_pos] = y_limit_pos
      
      p=ggplot(data=volvano_plot, aes(x=log2FoldChange, y=padj, color = color, label=lable_name)) + 
        geom_point() + 
        scale_color_manual(values=c("green", "grey", "red")) +
        theme(legend.position = "none") +
        geom_text_repel(aes(label = lable_name),size = 3, colour = 'black') 
      
      print(p)
    }
  }

  vsd <- varianceStabilizingTransformation(dds)  # transfer data using a variance stabilizing transformation

  ntd <- normTransform(dds)

  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:20]
  df <- as.data.frame(colData(dds)[,c("factor_1","factor_2")])
  colnames(df) <- c(factor_1,factor_2)
  pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
           cluster_cols=TRUE, annotation_col=df)
  
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:nrow(cts)]
  df <- as.data.frame(colData(dds)[,c("factor_1","factor_2")])
  colnames(df) <- c(factor_1,factor_2)
  pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
           cluster_cols=TRUE, annotation_col=df)
  
  
  sampleDists <- dist(t(assay(vsd)))
  
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$factor_1, vsd$factor_2, sep="_")     
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  

  pcaData <- plotPCA(vsd, intgroup=c("factor_2", "factor_1"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  ggplot(pcaData, aes(PC1, PC2, color=factor_2, shape=factor_1)) +
    geom_point(size=3) +
    labs(shape=factor_1, color=factor_2) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
  
  
} else {
  print("Only one or two factors are avaiable in the script")
}
  


