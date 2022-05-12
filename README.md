# EZ2RNAseq for RNA-seq data analysis
# The purpose of the pipeline is to let R packages easy to be used by researchers who do not have bioinformatics background. 
# DESeq2, limma, pasilla, pheatmap, vsn, RColorBrewer, ggplot2 and hexbin packages are used in the script.
# Before using the script, please install R and Rstudio in your computer. 
# 1. Prepare your reads data in ‘data.csv’ file and design information in ‘design.csv’ file. If you have only one factor in your experiment design, you can delete the column ‘C’ in ‘design.csv’. You can use Google Sheets or EXCEL to edit .csv files. If you want to annotate gene names in your volcano plot, you need to input gene names in ‘lable.csv’.
# 2. Open ‘EZ2RNAseq.R’ in Rstudio. Move cursor to command windows, press ‘Control+A’ on Windows or ‘command+A’ on Mac to choose all of the commands and click ‘run’.
