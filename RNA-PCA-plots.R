# for post esophagus samples -------------

# import counts table and and metadata file into R
# DeSeq2 uses non-normalized data

library(pasilla)

# start importing
setwd("~/d_abalone_caxc/RNAseq_analysis/diff_exp")

# import gene counts matrix
pasCts<- read.delim("Trinity_trans.gene.counts.annotated.matrix", row.names = 1)
cts <- as.matrix(pasCts,sep="\t")
# convert decimals to integers 
rounded_cts <- round(cts)

# import sample metadata
coldata <- read.delim("deseq_sample_data.txt", row.names = 1)
coldata$condition <- factor(coldata$condition)
coldata$tissue <- factor(coldata$tissue)

# check data
head(rounded_cts,2)
coldata

# create DeSeq2 object
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = rounded_cts,
                              colData = coldata,
                              design = ~ condition)
dds

# subset to just PE samples
PE_dds <- dds[, dds$tissue == "post esophagus"]
PE_dds

# pre-filter genes
smallestGroupSize <- 5
keep <- rowSums(counts(PE_dds) >= 10) >= smallestGroupSize
PE_dds <- PE_dds[keep,]

# run DESeq2
#by default, R will choose a reference level for factors based on alphabetical order
PE_DS <- DESeq(PE_dds)
# get the results
PE_res <- results(PE_DS)
PE_res

# filter to just significant results with log fold change greater than abs value of 2
# also omit NAs
PE_sigs <- na.omit(PE_res)
PE_sigs <- PE_sigs[PE_sigs$padj < 0.05,]
PE_sigs.df <- as.data.frame(PE_sigs)
PE_sigs.df <- PE_sigs.df[(abs(PE_sigs.df$log2FoldChange) > 2),]

# write table of significant transcript IDs
write.csv(PE_sigs.df, "PE_sigs_df.csv")

# make PCA plot of PE samples colored by condition
# variance stabilizing transformation (vst) -- remove the dependence of the variance on the mean
PE_vsd <- vst(PE_DS, blind=FALSE)
# first look
plotPCA(PE_vsd, intgroup=c("condition"))
# make more annotated plot
PE_pcaData <- plotPCA(PE_vsd, intgroup=c("condition"), returnData=TRUE)
PE_percentVar <- round(100 * attr(PE_pcaData, "percentVar"))
library(ggplot2)
ggplot(PE_pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",PE_percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",PE_percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_manual(values=c("darkcyan", "darkorange2")) +
  ggtitle("Post Esophagus Gene Expression PCA")


# for digestive gland samples ---------------

# import counts table and and metadata file into R
# DeSeq2 uses non-normalized data

library(pasilla)

# start importing
setwd("~/d_abalone_caxc/RNAseq_analysis/diff_exp")

# import gene counts matrix
pasCts<- read.delim("Trinity_trans.gene.counts.annotated.matrix", row.names = 1)
cts <- as.matrix(pasCts,sep="\t")
# convert decimals to integers 
rounded_cts <- round(cts)

# import sample metadata
coldata <- read.delim("deseq_sample_data.txt", row.names = 1)
coldata$condition <- factor(coldata$condition)
coldata$tissue <- factor(coldata$tissue)

# check data
head(rounded_cts,2)
coldata

# create DeSeq2 object
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = rounded_cts,
                              colData = coldata,
                              design = ~ condition)
dds

# subset to just PE samples
DG_dds <- dds[, dds$tissue == "digestive gland"]
DG_dds

# pre-filter genes
smallestGroupSize <- 5
keep <- rowSums(counts(DG_dds) >= 10) >= smallestGroupSize
DG_dds <- DG_dds[keep,]

# run DESeq2
#by default, R will choose a reference level for factors based on alphabetical order
DG_DS <- DESeq(DG_dds)
# get the results
DG_res <- results(DG_DS)
DG_res

# filter to just significant results with log fold change greater than abs value of 2
# also omit NAs
DG_sigs <- na.omit(DG_res)
DG_sigs <- DG_sigs[DG_sigs$padj < 0.05,]
DG_sigs.df <- as.data.frame(DG_sigs)
DG_sigs.df <- DG_sigs.df[(abs(DG_sigs.df$log2FoldChange) > 2),]

# write table of significant transcript IDs
write.csv(DG_sigs.df, "DG_sigs_df.csv")

# make PCA plot of DG samples colored by condition
# variance stabilizing transformation (vst) -- remove the dependence of the variance on the mean
DG_vsd <- vst(DG_DS, blind=FALSE)
# first look
plotPCA(DG_vsd, intgroup=c("condition"))
# make more annotated plot
DG_pcaData <- plotPCA(DG_vsd, intgroup=c("condition"), returnData=TRUE)
DG_percentVar <- round(100 * attr(DG_pcaData, "percentVar"))
library(ggplot2)
ggplot(DG_pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",PE_percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",PE_percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_manual(values=c("darkcyan", "darkorange2")) +
  ggtitle("Digestive Gland Gene Expression PCA")


# for comparing differentially expressed transcripts between tissue types ---------

# copy lists from PE_sigs_df.csv and DG_sigs_df.csv into Jvenn