# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
# Install writexl package
install.packages("writexl")

# Then load the library
library(writexl)
library(GEOquery)
library(Biobase)
library(GEOquery)
library(limma)
library(umap)
library(DESeq2)
# load series and platform data from GEO########### normalization methods: limma (common for both microarray and RNAseq) ############## However Deseq2 EdgeR
# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE245006", "file=GSE245006_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

# load gene annotations 
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID

# sample selection


##When you create a factor from numbers using gs <- factor(sml),
# R automatically sorts the unique values in ascending order to create the default levels. So:
# 0 would be the first level 0 → "Amsacrine" (first name in groups)
# 1 would be the second level 1 → "Control-DMSO" (second name)



# sample selection
gsms <- paste0("XXXXXXXXXXXXXXXXXXXXXXXXXX1XXXXXXXXXXXXXXXXXXXXXX0",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX1XXX1XXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXX1XX1XXXXXXXXXXXXXXXXXXXXXX0XXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX1XXX1XXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XX")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
tbl <- tbl[ ,sel]

# group membership for samples
gs <- factor(sml)
groups <- make.names(c("Amsacrine","Control-DMSO"))
levels(gs) <- groups
sample_info <- data.frame(Group = gs, row.names = colnames(tbl))

# pre-filter low count genes
# keep genes with at least N counts > 10, where N = size of smallest group
keep <- rowSums( tbl >= 10 ) >= min(table(gs))
tbl <- tbl[keep, ]

ds <- DESeqDataSetFromMatrix(countData=tbl, colData=sample_info, design= ~Group)

ds <- DESeq(ds, test="Wald", sfType="poscount")

# extract results for top genes table
r <- results (ds, contrast=c("Group", groups[1], groups[2]), alpha=0.05, pAdjustMethod ="fdr")

tT <- r[order(r$padj)[1:250],] 
tT <- merge(as.data.frame(tT), annot, by=0, sort=F)

tT <- subset(tT, select=c("GeneID","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean","Symbol","Description"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

plotDispEsts(ds, main="GSE245006 Dispersion Estimates")

# create histogram plot of p-values
hist(r$padj, breaks=seq(0, 1, length = 21), col = "grey", border = "white", 
     xlab = "", ylab = "", main = "GSE245006 Frequencies of padj-values")

# volcano plot
old.pal <- palette(c("#00BFFF", "#FF3030")) # low-hi colors
par(mar=c(4,4,2,1), cex.main=1.5)
plot(r$log2FoldChange, -log10(r$padj), main=paste(groups[1], "vs", groups[2]),
     xlab="log2FC", ylab="-log10(Padj)", pch=20, cex=0.5)
with(subset(r, padj<0.05 & abs(log2FoldChange) >= 0),
     points(log2FoldChange, -log10(padj), pch=20, col=(sign(log2FoldChange) + 3)/2, cex=1))
legend("bottomleft", title=paste("Padj<", 0.05, sep=""), legend=c("down", "up"), pch=20,col=1:2)

# MD plot
par(mar=c(4,4,2,1), cex.main=1.5)
plot(log10(r$baseMean), r$log2FoldChange, main=paste(groups[1], "vs", groups[2]),
     xlab="log10(mean of normalized counts)", ylab="log2FoldChange", pch=20, cex=0.5)
with(subset(r, padj<0.05 & abs(log2FoldChange) >= 0),
     points(log10(baseMean), log2FoldChange, pch=20, col=(sign(log2FoldChange) + 3)/2, cex=1))
legend("bottomleft", title=paste("Padj<", 0.05, sep=""), legend=c("down", "up"), pch=20,col=1:2)
abline(h=0)
palette(old.pal) # restore palette

################################################################
#   General expression data visualization
dat <- log10(counts(ds, normalized = T) + 1) # extract normalized counts

# box-and-whisker plot
lbl <- "log10(raw counts + 1)"
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
boxplot(dat[,ord], boxwex=0.6, notch=T, main="GSE245006", ylab="lg(norm.counts)", outline=F, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# UMAP plot (multi-dimensional scaling)
library(umap)
dat <- dat[!duplicated(dat), ] # first remove duplicates
par(mar=c(3,3,2,6), xpd=TRUE, cex.main=1.5)
ump <- umap(t(dat), n_neighbors = 4, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=4", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=groups, pch=20,
       col=1:length(groups), title="Group", pt.cex=1.5)

###############################################################################################################################################



# Get normalized counts
normalized_counts <- counts(ds, normalized=TRUE)

# Convert to data frame and add gene symbols
# First get the gene annotations for your normalized counts
normalized_counts_df <- as.data.frame(normalized_counts)
normalized_counts_df$GeneID <- rownames(normalized_counts_df)

# Merge with annotations to get gene symbols and descriptions
normalized_counts_with_info <- merge(normalized_counts_df, annot, by="GeneID", sort=FALSE)

# Write to Excel file
library(writexl)
write_xlsx(normalized_counts_with_info, "normalized_counts_with_annotations.xlsx")

# If you also want to save the differential expression results:
# Get the complete DESeq2 results
deseq_results <- as.data.frame(r)
deseq_results$GeneID <- rownames(deseq_results)

# Merge with annotations
deseq_results_with_info <- merge(deseq_results, annot, by="GeneID", sort=FALSE)

# Write DESeq2 results to Excel
write_xlsx(deseq_results_with_info, "deseq2_results_with_annotations.xlsx")


sample_info$sample_names <- rownames(sample_info)
write_xlsx(sample_info, "sample2.xlsx")
###################################################################################
# Load required package for heatmap
if (!require("pheatmap")) {
  install.packages("pheatmap")
}
library(pheatmap)

# Order results by log2FoldChange
ordered_results <- r[order(r$log2FoldChange), ]  # This will sort from smallest to largest log2FC
top_50_genes <- rownames(ordered_results)[1:50]

# Extract normalized counts for top 50 genes
norm_counts <- counts(ds, normalized=TRUE)
top_50_counts <- norm_counts[top_50_genes, ]

# Log transform the counts
log2_counts <- log2(top_50_counts + 1)

# Get gene symbols for the top 50 genes
top_50_info <- annot[top_50_genes, ]
row_labels <- top_50_info$Symbol

# Create annotation for columns
column_annotation <- data.frame(Group = gs)
rownames(column_annotation) <- colnames(log2_counts)

# First, let's check the actual levels in your group factor
print(levels(gs))

# Now modify the color scheme to match exactly
ann_colors <- list(Group = c("Control.DMSO" = "#FF3030", "Amsacrine" = "#00BFFF"))  # Note the dot instead of hyphen

# Create annotation for columns
column_annotation <- data.frame(Group = gs)
rownames(column_annotation) <- colnames(log2_counts)

# Create heatmap
pheatmap(log2_counts,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = column_annotation,
         annotation_colors = ann_colors,
         labels_row = row_labels,
         main = "Top 50 DEGs (Ordered by log2FoldChange)",
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         fontsize_row = 8,
         height = 10,
         width = 8)