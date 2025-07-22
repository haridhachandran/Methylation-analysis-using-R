# ===============================================================
# Title: Methylation Array Analysis using Illumina 450K Data
# Author: Haridha Sree C
# Updated: July 2025
# Description: Workflow for quality control and preprocessing of
#              Illumina 450k methylation array data using minfi.
# ===============================================================

# --- Step 1: Install Bioconductor & Required Packages ---
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")

# Install methylation analysis package
BiocManager::install("methylationArrayAnalysis")

# --- Step 2: Load Required Libraries ---
library(knitr)
library(limma)
library(minfi)
library(methylationArrayAnalysis)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)

# --- Step 3: Load Example Data ---
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
list.files(dataDirectory, recursive = TRUE)

# --- Step 4: Read Sample Sheet & Raw Data (.idat files) ---
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheet.csv")
rgSet <- read.metharray.exp(targets = targets)

# Assign readable sample names
targets$ID <- paste(targets$Sample_Group, targets$Sample_Name, sep = ",")
sampleNames(rgSet) <- targets$ID

# View RGChannelSet
print(rgSet)

# --- Step 5: Quality Control - Detection P-values ---
detP <- detectionP(rgSet)
head(detP)

# Bar plot of mean detection p-values
pal <- brewer.pal(8, "Dark2")
par(mfrow = c(1, 2))

# Default scale
barplot(colMeans(detP), col = pal[factor(targets$Sample_Group)], las = 2,
        cex.names = 0.8, ylab = "Mean Detection p-values")
abline(h = 0.05, col = "red")
legend("topleft", legend = levels(factor(targets$Sample_Group)), fill = pal, bg = "white")

# Zoomed in (threshold scale)
barplot(colMeans(detP), col = pal[factor(targets$Sample_Group)], las = 2,
        cex.names = 0.8, ylim = c(0, 0.002), ylab = "Mean Detection p-values")
abline(h = 0.05, col = "red")
legend("topleft", legend = levels(factor(targets$Sample_Group)), fill = pal, bg = "white")

# --- Step 6: Save Quality Control Report ---
qcReport(rgSet,
         sampNames = targets$ID,
         sampGroups = targets$Sample_Group,
         pdf = "qcReport.pdf")

# --- Step 7: Sample Filtering based on Quality ---
# Remove samples with mean detection p-values > 0.05
keep_samples <- colMeans(detP) < 0.05
rgSet <- rgSet[, keep_samples]
print(rgSet)
