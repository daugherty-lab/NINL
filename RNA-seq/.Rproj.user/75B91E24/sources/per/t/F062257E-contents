
# R script for DESeq2 processing of NIN-KO and NIN-KO_induced
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd('..')
#### installing and loading packages ####

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos='http://cran.us.r-project.org')

knitr::opts_chunk$set(tidy=FALSE, cache=TRUE,
                      dev="png",
                      message=FALSE, error=FALSE, warning=TRUE)
options(scipen = 5) # scientific notation for plotting later

packages <- c("reshape", "ggplot2", "ggrepel", "RColorBrewer", "pheatmap", "data.table")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())), repos='http://cran.us.r-project.org')  
}
install.packages("data.table", dependencies=TRUE)
BiocManager::install("DESeq2")
BiocManager::install("DEGreport")
BiocManager::install("GenomicFeatures")

library(reshape)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
suppressPackageStartupMessages(library(DESeq2))
library(pheatmap)
suppressPackageStartupMessages(library('GenomicFeatures'))
library(tximport)
library(readr)
library(tximportData)

##  ----config---- 
species <- "NIN-KO_induced_to_NIN-KO"

## ---- making an tx2gene for hg38 ------------ 
txdb <-makeTxDbFromGFF("/Users/bryantcao/Desktop/Homo_sapiens.GRCh38.100.gtf")
keytypes(txdb)
k <- keys(txdb, keytype = "TXNAME")
df <- select(txdb, keys = k,  columns = "GENEID", keytype = "TXNAME")
tx2gene <- df[, 1:2]
head(tx2gene)
write.csv(tx2gene, paste0("out/tx/", species, "_tx2gene.csv"), row.names = FALSE)

## ----txiSetup------------------------------------------------------------
dir <- system.file("/chrisRNA/quants/", package="tximportData")
dir
#samples <- read.table(file.path(dir, paste0("Volumes/MiniDrive/Dropbox/daugherty-lab/NINL/lists/", species, ".txt")), header=TRUE)
samples <- read.table(file.path(dir, paste0("Volumes/MiniDrive/Dropbox/daugherty-lab/NINL/lists/", "all", ".txt")), header=TRUE)
samples$Condition
rownames(samples) <- samples$Run
samples
## ----txiFiles------------------------------------------------------------
files <- file.path(paste0("/Volumes/MiniDrive/Dropbox/daugherty-lab/NINL/chrisRNA/quants/", samples$Run, ".sf"))
names(files) <- samples$run
tx2gene <- read_csv(file.path(paste0("/Volumes/MiniDrive/Dropbox/daugherty-lab/NINL/out/tx/", species, "_tx2gene.csv")))

## ----tximport, results="hide"--------------------------------------------
txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = TRUE)

## ----txi2dds, results="hide"---------------------------------------------

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ Condition)
dds <- ddsTxi

## ----prefilter-----------------------------------------------------------
#keep <- rowSums(counts(dds)) >= 1
#dds <- dds[keep,]

## ----relevel-------------------------------------------------------------
dds$Condition <- relevel(dds$Condition, ref = "NIN-KO")

## ----droplevels----------------------------------------------------------
#dds$condition <- droplevels(dds$condition)

## ----deseq---------------------------------------------------------------
dds <- DESeq(dds)

res <- results(dds)

# PICK ONLY ONE SET OF CONDITIONS TO RUN AT A TIME
res <- results(dds, contrast=c("Condition","NIN-KO_induced","NIN-KO"))

res


## ----export, eval=FALSE--------------------------------------------------
# RUN THE APPROPRIATE ONE##########################

write.csv(as.data.frame(res), file=("out/processed/Condition_NIN-KO_induced_vs_NIN-KO.csv"))

#----volcano plot - expression threshold focused ----

## Obtain logical vector regarding whether padj values are less than 0.1
threshold_OE <- res$padj < 1 
threshold_OE
## Determine the number of TRUE values
length(which(threshold_OE))

res$threshold <- threshold_OE 

df_res = as.data.frame(res)

df_res

## Volcano plot
ggplot(df_res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  ggtitle(species) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,100)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  

______________________

vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)
# 
# rld <- rlog(dds, blind = FALSE)
# head(assay(rld), 3)
# 
# library("dplyr")
# library("ggplot2")
# 
# dds <- estimateSizeFactors(dds)
# 
# df <- bind_rows(
#   as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
#     mutate(transformation = "log2(x + 1)"),
#   as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
#   as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
# 
# colnames(df)[1:2] <- c("x", "y")  
# 
# ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
#   coord_fixed() + facet_grid( . ~ transformation)  
sampleDists <- dist(t(assay(vsd)))
sampleDists

library("pheatmap")
library("RColorBrewer")  

nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$Run )
colnames(sampleDistMatrix) <- paste( vsd$Run )
colors <- colorRampPalette( rev(brewer.pal(12, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
sampleDistMatrix

plotPCA(vsd, intgroup = c("Condition", "Run"))


####----not working due to too many shapes----
 pcaData <- plotPCA(vsd, intgroup = c( "Condition", "Run"), returnData = TRUE)
 pcaData  
 percentVar <- round(100 * attr(pcaData, "percentVar"))
 p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Run, shape = Condition)) +
   scale_fill_manual(values = mycolors) +
   geom_point(size =5, alpha = 0.6) +
   xlab(paste0("PC1: ", percentVar[1], "% variance")) +
   ylab(paste0("PC2: ", percentVar[2], "% variance")) +
   coord_fixed()

ggsave('test_PCA.png', width = 35, height = 20, units = 'cm')

