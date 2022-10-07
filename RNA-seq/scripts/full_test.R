

# R script for DESeq2 processing of 1

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
species <- 'hg38'

## ---- making an tx2gene for hg38 ------------ 
txdb <-makeTxDbFromGFF("Homo_sapiens.GRCh38.100.gff")
keytypes(txdb)
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k,  columns = "TXNAME", keytype = "GENEID")
tx2gene <- df[, 2:1]
head(tx2gene)
write.csv(tx2gene, paste0("out/tx/", species, "_tx2gene.csv"), row.names = FALSE)

## ----txiSetup------------------------------------------------------------
#library("tximport")
#library("readr")
library("tximportData")
dir <- system.file("/chrisRNA/quants/", package="tximportData")
samples <- read.table(file.path(dir, paste0("Volumes/MiniDrive/Dropbox/daugherty-lab/Helitrons/DESeq2/lists/", species, ".txt")), header=TRUE)
samples$Condition
rownames(samples) <- samples$Run
samples
## ----txiFiles------------------------------------------------------------
files <- file.path(paste0("/Volumes/MiniDrive/Dropbox/daugherty-lab/Helitrons/DESeq2/quants/", samples$Run, ".sf"))
files
names(files) <- samples$run
tx2gene <- read_csv(file.path(paste0("/Volumes/MiniDrive/Dropbox/daugherty-lab/Helitrons/DESeq2/tx/", species, "_tx2gene.csv")))

## ----tximport, results="hide"--------------------------------------------
txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreAfterBar = TRUE)

## ----txi2dds, results="hide"---------------------------------------------
library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ Condition)
dds <- ddsTxi


