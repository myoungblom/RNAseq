require(DESeq2)
require(dplyr)
require(ggplot2)
require(vsn)
require(RColorBrewer)
require(pheatmap)
require(apeglm)
require(genefilter)
require(hexbin)

#####
# Analysis of RNAseq data using DESeq2 from HTSeq count files.
# Requires: metadatafile, HTSeq count files 
#####

# get metadata file from working directory
setwd("/Users/madi/Desktop/06.18.2019_MT345RNAseq")
metadataFile <- grep("metadata.txt", list.files("."), value = TRUE)
sampleData <- read.delim(metadataFile, sep="", header = TRUE)

# reorganize metadata file for input to DESeq2
sampleTable <- data.frame(SampleName = sampleData$LibraryName,
                          CountsFile = sampleData$CountFile,
                          Genotype = factor(sampleData$Genotype),
                          Condition = factor(sampleData$Condition),
                          SampleID = sampleData$SampleID)

# create sample table for DESeq2 input
DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = ".",
                                          design = ~ Genotype+Condition)

# add annotation column to DESeq2Table from HTSeq counts file 
# that contains gene names in the second column
countsAnnot <- read.delim("/Users/madi/Desktop/06.18.2019_MT345RNAseq/33-A-P_HTSeqCounts_geneNames.txt", 
                         sep="", header=FALSE, col.names = c("gene_id","gene_name"))
mcols(DESeq2Table) <- cbind(mcols(DESeq2Table), countsAnnot$gene_name) # not sure if this is working properly

# remove rows with 0 counts
DESeq2Table <- DESeq2Table[rowSums(counts(DESeq2Table)) > 1,]

# quality control: heatmap
rld <- rlogTransformation(DESeq2Table, blind=FALSE) # log transformation
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists) # create distance matrix for heatmap
rownames(sampleDistMatrix) <- paste(rld$SampleID)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(100)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, legend = TRUE)

# quality control: PCA
DESeq2::plotPCA(rld, intgroup=c("SampleID"))

# differential expression analysis
DESeq2Table <- DESeq(DESeq2Table)
result <- results(DESeq2Table, alpha = 0.05, lfcThreshold = 1)
AvsEresult <- results(DESeq2Table, contrast =c("Genotype","Ancestral","Evolved"), alpha = 0.05, lfcThreshold = 1)

# num of genes from pellicle vs biofilm comparison with p value < 0.05
sum(result$padj <0.05, na.rm = TRUE)
sum(AvsEresult$pvalue <0.05, na.rm=TRUE)

# subset and sort to get genes with strongest up/down regulation in PvsB
subSet <- subset(result, padj < 0.05)
topDown <- head(subSet[ order(subSet$log2FoldChange), ], n=10)
topUp <- head(subSet[ order(subSet$log2FoldChange, decreasing = TRUE), ], n=10)

# MA plot
resultsNames(DESeq2Table)
maData <- lfcShrink(DESeq2Table, coef = "Condition_Ancestral_vs_Evolved", type="apeglm")
plotMA(maData, ylim = c(-4, 12))

# gene clustering heat map of all genes
allGenes <- (order(rowVars(assay(rld)), decreasing = TRUE))
mat <- assay(rld)[ allGenes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, c("Genotype", "Condition")])
pheatmap(mat, show_rownames = FALSE, annotation_col = anno)

# gene clustering heat map of most variable genes
topVarGenes <- head((order(rowVars(assay(rld)), decreasing = TRUE)), 40)
mat2 <- assay(rld)[ topVarGenes, ]
mat2 <- mat2 - rowMeans(mat2)
pheatmap(mat2, show_rownames = FALSE, annotation_col = anno)
