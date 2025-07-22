# Methylation-analysis-using-R
DNA Methylation Array QC and Preprocessing using Illumina 450k Arrays in R


**👩‍🔬 Overview:**
This project implements a standardized pipeline to perform loading, quality control, and preliminary filtering of Illumina HumanMethylation450k array datasets using the minfi and methylationArrayAnalysis packages in R.


**The analysis flow:**

Loads .idat files via a sample sheet

Performs detection p-value-based quality assessment

Filters low-quality samples

Generates visual diagnostics and a PDF quality control report

The approach utilizes public sample datasets available within the methylationArrayAnalysis Bioconductor package and serves as a template for user-generated data analysis.



**🎯 Key Features:**

📥 Automatic data loading from .idat arrays via CSV SampleSheet

🔍 Basic quality control using detection P-values

📊 Group-wise summary visualization with color-coded plots

📄 Exportable PDF quality control report

🔧 Preprocessing steps for downstream DMR/β-value analysis


**🛠 Technology Stack:**
Language: R

Packages: minfi, limma, methylationArrayAnalysis, missMethyl, DMRcate, illuminaHumanMethylation450kanno, Gviz


**📂 Input Requirements:**

A properly formatted SampleSheet.csv

Corresponding .idat raw files

A known experimental design (e.g., case-control groups)


**📤 Output:**

Mean detection p-value bar plots before and after sample removal

qcReport.pdf: Detailed quality report

RGChannelSet object with high-quality samples only
