
# R script for DESeq2 processing of WT and NINL-KO
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
species <- "NINL-KO_to_WT"

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
samples <- read.table(file.path(dir, paste0("Volumes/MiniDrive/Dropbox/daugherty-lab/NINL/lists/", species, ".txt")), header=TRUE)
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
dds$Condition <- relevel(dds$Condition, ref = "WT")

## ----droplevels----------------------------------------------------------
#dds$condition <- droplevels(dds$condition)

## ----deseq---------------------------------------------------------------
dds <- DESeq(dds)

res <- results(dds)

# PICK ONLY ONE SET OF CONDITIONS TO RUN AT A TIME
res <- results(dds, contrast=c("Condition","NINL-KO","WT"))

res


## ----export, eval=FALSE--------------------------------------------------
# RUN THE APPROPRIATE ONE##########################

write.csv(as.data.frame(res), file=("out/processed/Condition_NINL-KO_vs_WT.csv"))

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
