##### install packages #####
ip <- as.data.frame(installed.packages())
ip <- ip$Package

if (sum(ip == "rstudioapi") == 0) {
  install.packages("rstudioapi")
}

if (sum(ip == "vsn") == 0) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("vsn")
}

if (sum(ip == "DESeq2") == 0) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("DESeq2")
}

if (sum(ip == "hexbin") == 0) {
  install.packages("hexbin")
}

if (sum(ip == "pheatmap") == 0) {
  install.packages("pheatmap")
}

if (sum(ip == "RColorBrewer") == 0) {
  install.packages("RColorBrewer")
}

if (sum(ip == "ggplot2") == 0) {
  install.packages("ggplot2")
}

if (sum(ip == "ggrepel") == 0) {
  install.packages("ggrepel")
}

###################################################################
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

if (sum(colnames(cts) %in% row.names(design)) == 0) {
  cts = as.data.frame(t(cts))
}

if (sum(colnames(cts) %in% row.names(design)) == 0) {
  print('Sample_ID in the metadata is not the same as that in the reads table')
}

keep = colnames(cts) %in% row.names(design)
cts = cts[,keep]
cts = cts[row.names(design)]

x = min(apply(cts,2,min))
if (x == 0 ){
  cts = cts+1
}

path_input = dirname(rstudioapi::getSourceEditorContext()$path)
path_output = path_input
file_list = list.files(path = ".")
path_output = paste0(path_input,'/results')
path_output_1 = path_output
n=1
name_1 = 'results'
while (sum(!(name_1 %in% file_list)) == 0) {
  path_output_1 = paste0(path_output,n)
  name_1 = c(name_1,paste0('results',n))
  n = n+1
}
dir.create(path_output_1)
setwd(path_output_1)


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
  
  i=2
    res2 <- results(dds, name=res[i])
    write.csv(res2,file = paste(res1[i],".csv"))
    
    # volcano plot
    lable_list = file.exists(paste0(path_input,'/lable.csv'))
    if (lable_list) {
      lable_list = file.info(paste0(path_input,'/lable.csv'))
      lable_list = lable_list$size
      
      if (lable_list != 0) {
        lable_list = as.matrix(read.csv(paste0(path_input,'/lable.csv'), header = F))
      } else {
        lable_list = NA
      }
      
    } else {
      lable_list = NA
    }
    
    
      volcano_plot <- as.data.frame(res2@listData)
      
      volcano_plot$gene_id = row.names(cts)
      volcano_plot$lable_name = NA
      volcano_plot$lable_name[volcano_plot$gene_id %in% lable_list] <- volcano_plot$gene_id[volcano_plot$gene_id %in% lable_list]
      
      volcano_plot$color <- 'grey'
      volcano_plot$color[volcano_plot$log2FoldChange >= 1 & volcano_plot$padj <= 0.05] = 'green'
      volcano_plot$color[volcano_plot$log2FoldChange <= -1 & volcano_plot$padj <= 0.05] = 'red'
      
      
      volcano_plot$padj <- -log10(volcano_plot$padj)
      volcano_plot = volcano_plot[!is.na(volcano_plot$padj) & !is.na(volcano_plot$padj),]
      
      if (is.na(lable_list)[1]) {
        ggplot(data=volcano_plot, aes(x=log2FoldChange, y=padj)) + 
          geom_point(aes(color = color)) + 
          scale_color_manual(values=c("green" = "green", "grey" = "grey", "red" = "red")) +
          theme(legend.position = "none") 
        
        ggsave(paste0(col_name[i-1],"_volcano.pdf"), width=4, height=4)
        
        x_limit_neg <- min(quantile(volcano_plot$log2FoldChange,.01, na.rm =T),-2)
        x_limit_pos <- max(quantile(volcano_plot$log2FoldChange,.99, na.rm =T),2)
        volcano_plot$log2FoldChange[volcano_plot$log2FoldChange < x_limit_neg] = x_limit_neg
        volcano_plot$log2FoldChange[volcano_plot$log2FoldChange > x_limit_pos] = x_limit_pos
        
        y_limit_pos <- max(quantile(volcano_plot$padj,.99, na.rm =T),2)
        volcano_plot$padj[volcano_plot$padj > y_limit_pos] = y_limit_pos
        
        ggplot(data=volcano_plot, aes(x=log2FoldChange, y=padj)) + 
          geom_point(aes(color = color)) + 
          scale_color_manual(values=c("green" = "green", "grey" = "grey", "red" = "red")) +
          theme(legend.position = "none") 
        
        ggsave(paste0(col_name[i-1],"_volcano_2.pdf"), width=4, height=4)
        
      } else {
        ggplot(data=volcano_plot, aes(x=log2FoldChange, y=padj, label=lable_name)) + 
          geom_point(aes(color = color)) + 
          scale_color_manual(values=c("green" = "green", "grey" = "grey", "red" = "red")) +
          theme(legend.position = "none") +
          geom_text_repel(aes(label = lable_name),size = 3, colour = 'black') 
        
        ggsave(paste0(col_name[i-1],"_volcano.pdf"), width=4, height=4)
        
        x_limit_neg <- min(quantile(volcano_plot$log2FoldChange,.01, na.rm =T),-2)
        x_limit_pos <- max(quantile(volcano_plot$log2FoldChange,.99, na.rm =T),2)
        volcano_plot$log2FoldChange[volcano_plot$log2FoldChange < x_limit_neg] = x_limit_neg
        volcano_plot$log2FoldChange[volcano_plot$log2FoldChange > x_limit_pos] = x_limit_pos
        
        y_limit_pos <- max(quantile(volcano_plot$padj,.99, na.rm =T),2)
        volcano_plot$padj[volcano_plot$padj > y_limit_pos] = y_limit_pos
        
        ggplot(data=volcano_plot, aes(x=log2FoldChange, y=padj, label=lable_name)) + 
          geom_point(aes(color = color)) + 
          scale_color_manual(values=c("green" = "green", "grey" = "grey", "red" = "red")) +
          theme(legend.position = "none") +
          geom_text_repel(aes(label = lable_name),size = 3, colour = 'black') 
        
        ggsave(paste0(col_name[i-1],"_volcano_2.pdf"), width=4, height=4)
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
  pdf("pheatmap_1.pdf", width=6, height=6)
  print(pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
                 cluster_cols=TRUE, annotation_col=df))
  
  dev.off()
  
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:nrow(cts)]
  df <- as.data.frame(colData(dds)[1])
  colnames(df) <- c(factor_1)
  
  pdf("pheatmap_2.pdf", width=6, height=6)
  print(pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=F,show_colnames=F,
                 cluster_cols=TRUE, annotation_col=df))
  dev.off()  
  
  ############### Heatmap of the sample-to-sample distances ###############
  sampleDists <- dist(t(assay(vsd)))
  
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$factor_1)     
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

  pdf("pheatmap_3.pdf", width=6, height=6)
  print(pheatmap(sampleDistMatrix,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,
                 col=colors))
  dev.off()
  
  
  ############### pincipal component plot of the samples ###############
  pcaData <- plotPCA(vsd, intgroup=c("factor_1"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  #pcaData <- avgdist(cts, sample = )
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  ggplot(pcaData, aes(PC1, PC2, color=factor_1)) +
    geom_point(size=2) +
    labs(shape=factor_1) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
  ggsave(("PCA.pdf"), width=6, height=6)
  
  ggplot(pcaData, aes(PC1, PC2, color=factor_1)) +
    geom_point(size=2) +stat_ellipse(type = "norm")+
    labs(shape=factor_1) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
  ggsave(("PCA_2.pdf"), width=6, height=6)
  
  
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
  res2 <- results(dds, name=res[4])
  write.csv(res2,file = paste(res1[4],".csv"))
  
  for (i in 2:3)   {
    res2 <- results(dds, name=res[i])
    write.csv(res2,file = paste(res1[i],".csv"))
    
    # volcano plot
    lable_list = file.exists(paste0(path_input,'/lable.csv'))
    if (lable_list) {
      lable_list = file.info(paste0(path_input,'/lable.csv'))
      lable_list = lable_list$size
      
      if (lable_list != 0) {
        lable_list = as.matrix(read.csv(paste0(path_input,'/lable.csv'), header = F))
      } else {
        lable_list = NA
      }
      
    } else {
      lable_list = NA
    }
    

      volcano_plot <- as.data.frame(res2@listData)
      
      volcano_plot$gene_id = row.names(cts)
      volcano_plot$lable_name = NA
      volcano_plot$lable_name[volcano_plot$gene_id %in% lable_list] <- volcano_plot$gene_id[volcano_plot$gene_id %in% lable_list]
      
      volcano_plot$color <- 'grey'
      volcano_plot$color[volcano_plot$log2FoldChange >= 1 & volcano_plot$padj <= 0.05] = 'green'
      volcano_plot$color[volcano_plot$log2FoldChange <= -1 & volcano_plot$padj <= 0.05] = 'red'
      
      
      volcano_plot$padj <- -log10(volcano_plot$padj)
      
      ggplot(data=volcano_plot, aes(x=log2FoldChange, y=padj, label=lable_name)) + 
        geom_point(aes(color = color)) + 
        scale_color_manual(values=c("green" = "green", "grey" = "grey", "red" = "red")) +
        theme(legend.position = "none") +
        geom_text_repel(aes(label = lable_name),size = 3, colour = 'black') 
      
      ggsave(paste0(col_name[i-1],"_volcano.pdf"), width=4, height=4)
      
      x_limit_neg <- min(quantile(volcano_plot$log2FoldChange,.01, na.rm =T),-2)
      x_limit_pos <- max(quantile(volcano_plot$log2FoldChange,.99, na.rm =T),2)
      volcano_plot$log2FoldChange[volcano_plot$log2FoldChange < x_limit_neg] = x_limit_neg
      volcano_plot$log2FoldChange[volcano_plot$log2FoldChange > x_limit_pos] = x_limit_pos
      
      y_limit_pos <- max(quantile(volcano_plot$padj,.99, na.rm =T),2)
      volcano_plot$padj[volcano_plot$padj > y_limit_pos] = y_limit_pos
      
      ggplot(data=volcano_plot, aes(x=log2FoldChange, y=padj, label=lable_name)) + 
        geom_point(aes(color = color)) + 
        scale_color_manual(values=c("green" = "green", "grey" = "grey", "red" = "red")) +
        theme(legend.position = "none") +
        geom_text_repel(aes(label = lable_name),size = 3, colour = 'black') 
      
      ggsave(paste0(col_name[i-1],"_volcano_2.pdf"), width=4, height=4)
    }
  

  vsd <- varianceStabilizingTransformation(dds)  # transfer data using a variance stabilizing transformation

  ntd <- normTransform(dds)

  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:20]
  df <- as.data.frame(colData(dds)[,c("factor_1","factor_2")])
  colnames(df) <- c(factor_1,factor_2)
  
  pdf("pheatmap_1.pdf", width=6, height=6)
  print(pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
                   cluster_cols=TRUE, annotation_col=df))

  dev.off()
  
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:nrow(cts)]
  df <- as.data.frame(colData(dds)[,c("factor_1","factor_2")])
  colnames(df) <- c(factor_1,factor_2)
  
  pdf("pheatmap_2.pdf", width=6, height=6)
  print(pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
                 cluster_cols=TRUE, annotation_col=df))
  dev.off()
  
  sampleDists <- dist(t(assay(vsd)))
  
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$factor_1, vsd$factor_2, sep="_")     
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  pdf("pheatmap_3.pdf", width=6, height=6)
  print(pheatmap(sampleDistMatrix,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,
                 col=colors))
  dev.off()

  pcaData <- plotPCA(vsd, intgroup=c("factor_2", "factor_1"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  ggplot(pcaData, aes(PC1, PC2, color=factor_2, shape=factor_1)) +
    geom_point(size=3) +
    labs(shape=factor_1, color=factor_2) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
  ggsave(("PCA.pdf"), width=4, height=4)
  
  ggplot(pcaData, aes(PC1, PC2, color=factor_2, shape=factor_1)) +
    geom_point(size=3) +
    labs(shape=factor_1, color=factor_2) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()+stat_ellipse(type = "t")
  ggsave(("PCA_2.pdf"), width=4, height=4)
  
} else {
  print("Only one or two factors are avaiable in the script")
}
  


