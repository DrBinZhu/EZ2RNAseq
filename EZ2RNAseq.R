###################################################################
# This version of DESeq2 script is edited by Dr. Bin Zhu at Virginia Commonwealth University
# Email: binzhu0824@gmail.com
###################################################################
###################################################################
##Start##

#load DESeq package;
library("DESeq2")

# set files path
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# Read in your counts table;
cts= read.csv("data.csv", header=TRUE, row.names=1)
design = read.csv("design.csv", header=TRUE, row.names=1)

opo <- dim(design)[2]

if (opo == 1) {
        
# change column names
col_name <- colnames(design)
factor_1 = col_name[1]

colnames(design) <- c("factor_1")

levels(design$factor_1)                                       

#   pick up a subgroup of samples   
#Sample <- design$Treatment == 'No_selenium'| design$Treatment == 'Selenium_from_B._subtilis'
#cts <- cts[,Sample]
#design <- design[Sample,]

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = design,
                              design= ~ factor_1)  
ref1 = toString(design$factor_1[1])

dds$factor_1 <- relevel(dds$factor_1, ref = ref1)   

#Pre-filtering
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]

dds <- DESeq(dds)


res=resultsNames(dds) # lists the coefficients

res1 <- gsub("factor_1", factor_1, res)
res1

n=length(res)

for (i in 2:n)   {
  res2 <- results(dds, name=res[i])
  write.csv(res2,file = paste(res1[i],".csv"))}




# Log fold change shrinkage for visualization and ranking
#resLFC <- lfcShrink(dds, coef="Concentration_Concentration_50_vs_Concentration_0", type="apeglm")
# Exploring and exporting results
# MA-plot
#plotMA(res1, ylim=c(-2,2))
#plotMA(resLFC, ylim=c(-2,2))




###              Data transformations and visualization         ###
vsd <- varianceStabilizingTransformation(dds)  # transfer data using a variance stabilizing transformation

# Effects of transformations on the variance   
# The figure below plots the standard deviation of the transformed data, across samples, against the mean, using the shifted logarithm transformation, the regularized log transformation and the variance stabilizing transformation. The shifted logarithm has elevated standard deviation in the lower count range, and the regularized log to a lesser extent, while for the variance stabilized data the standard deviation is roughly constant along the whole dynamic range.

# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
library("hexbin")
meanSdPlot(assay(ntd))


# Heatmap of the count matrix           ###
library("pheatmap")  # using R script 'hclust' to test distance 
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[1])
colnames(df) <- c(factor_1)
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)





# Heatmap of the sample-to-sample distances

sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$factor_1)     
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


# pincipal component plot of the samples
library("ggplot2")
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

} else if (opo == 2) {
  # change column names
  col_name <- colnames(design)
  factor_1 = col_name[1]
  factor_2 = col_name[2]
  
  colnames(design) <- c("factor_1","factor_2")
  
  levels(design$factor_1)                                       
  levels(design$factor_2)                                           
  
  #   pick up a subgroup of samples   
  #Sample <- design$Treatment == 'No_selenium'| design$Treatment == 'Selenium_from_B._subtilis'
  #cts <- cts[,Sample]
  #design <- design[Sample,]
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = design,
                                design= ~ factor_1 + factor_2 + factor_2:factor_1 )  
  ref1 = toString(design$factor_1[1])
  ref2 = toString(design$factor_2[1])
  
  dds$factor_1 <- relevel(dds$factor_1, ref = ref1)   
  dds$factor_2 <- relevel(dds$factor_2, ref = ref2)   
  
  #Pre-filtering
  #keep <- rowSums(counts(dds)) >= 10
  #dds <- dds[keep,]
  
  dds <- DESeq(dds)
  
  
  res=resultsNames(dds) # lists the coefficients
  
  res1 <- gsub("factor_1", factor_1, res)
  res1 <- gsub("factor_2", factor_2, res1)
  res1
  
  n=length(res)
  
  for (i in 2:n)   {
    res2 <- results(dds, name=res[i])
    write.csv(res2,file = paste(res1[i],".csv"))}
  
  
  
  
  # Log fold change shrinkage for visualization and ranking
  #resLFC <- lfcShrink(dds, coef="Concentration_Concentration_50_vs_Concentration_0", type="apeglm")
  # Exploring and exporting results
  # MA-plot
  #plotMA(res1, ylim=c(-2,2))
  #plotMA(resLFC, ylim=c(-2,2))
  
  
  
  
  ###              Data transformations and visualization         ###
  vsd <- varianceStabilizingTransformation(dds)  # transfer data using a variance stabilizing transformation
  
  # Effects of transformations on the variance   
  # The figure below plots the standard deviation of the transformed data, across samples, against the mean, using the shifted logarithm transformation, the regularized log transformation and the variance stabilizing transformation. The shifted logarithm has elevated standard deviation in the lower count range, and the regularized log to a lesser extent, while for the variance stabilized data the standard deviation is roughly constant along the whole dynamic range.
  
  # this gives log2(n + 1)
  ntd <- normTransform(dds)
  library("vsn")
  library("hexbin")
  meanSdPlot(assay(ntd))
  
  
  # Heatmap of the count matrix           ###
  library("pheatmap")  # using R script 'hclust' to test distance 
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:20]
  df <- as.data.frame(colData(dds)[,c("factor_2","factor_1")])
  colnames(df) <- c(factor_1,factor_2)
  pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
           cluster_cols=TRUE, annotation_col=df)
  
  
  
  
  
  # Heatmap of the sample-to-sample distances
  
  sampleDists <- dist(t(assay(vsd)))
  
  library("RColorBrewer")
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$factor_1, vsd$factor_2, sep="_")     
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  
  
  # pincipal component plot of the samples
  library("ggplot2")
  pcaData <- plotPCA(vsd, intgroup=c("factor_2", "factor_1"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  #pcaData <- avgdist(cts, sample = )
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
  


