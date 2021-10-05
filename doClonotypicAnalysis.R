#!/usr/bin/Rscript --vanilla

# Aim:
# summarize clonotype data, plot data.


# Usage: 
# Rscript doClonotypicAnalysis.R


library(tidyverse)
source("./R/functions.R")


# inputs
meta_path <- "./zebrafish/meta/samples_zf_RNAseq_Mariana.xlsx"
sheet_name <- "main_editCP_missing_sample"
results_path <- "./zebrafish/results/"
sample_name_column_in_meta <- "sample_name"
sample_name_column_in_data <- "sample"

arrangement_type <- "productive"
receptor_type <- "TRB"


file <- paste0(results_path, "/clonotypes_raw_wMeta_wRatios.TRB.productive.tsv")

clonotypes <- data.table::fread(file)

# Optional: add extra column with custum names to metadata
metadata <- readxl::read_xlsx(path = meta_path, sheet = sheet_name) %>% 
  distinct()

# =================================================================================================
# Plot the clonotypes 
# =================================================================================================



# for this particular analysis we are interested in plotting subgroups of the data: 
# 1. Get summarized data: 
uniqueRearr <- clonotypes %>% 
  totalUniqueRearrPerSamplePerReceptor(
    receptor_column = "allVHitsWithScore", 
    receptor_type = "TRB", 
    sample_column = "sample", 
    outdir = results_path
  )

# 1.1 plot number of unique clonotypes
# 1.1.2 Filter (in) the MYC samples: 
# 1.1.2.1 vertical 
uniqueRearr %>% 
  filter(str_detect(`sample`, "Myc")) %>% 
  plotBoxplotDotJitter(
    grouping = "Manuscript notation", 
    outname = "TransientLines", 
    ylab  = paste0("#clonotypes(",receptor_type,")"),
    subDir = arrangement_type, 
    outdir = results_path,
    meta_path = meta_path,
    sheet_name = sheet_name
  )
# 1.1.2.2 horizontal
uniqueRearr %>% 
  filter(str_detect(`sample`, "Myc")) %>% 
  plotBoxplotDotJitterHz(
    grouping = "Manuscript notation", 
    outname = "TransientLines",
    ylab  = paste0("#clonotypes(",receptor_type,")"),
    subDir = arrangement_type, 
    outdir = results_path,
    meta_path = meta_path,
    sheet_name = sheet_name
  )

# 1.1.3 Filter (in) the WKM_vs_stable
uniqueRearr  %>% 
  filter(!str_detect(`sample`, "Myc")) %>% 
  filter(!str_detect(`sample`, "Thy")) %>% 
  filter(!str_detect(`sample`, "WKM-CG1-P2mut")) %>% 
  plotBoxplotDotJitter(
    grouping = "Manuscript notation", 
    outname = "WKM_vs_stable", 
    ylab  = paste0("#clonotypes(",receptor_type,")"),
    subDir = "productive",
    outdir = results_path,
    meta_path = meta_path,
    sheet_name = sheet_name
  )

uniqueRearr  %>% 
  filter(!str_detect(`sample`, "Myc")) %>% 
  filter(!str_detect(`sample`, "Thy")) %>% 
  filter(!str_detect(`sample`, "WKM-CG1-P2mut")) %>% 
  plotBoxplotDotJitterHz(
    grouping = "Manuscript notation", 
    outname = "WKM_vs_stable", 
    ylab  = paste0("#clonotypes(",receptor_type,")"),
    subDir = "productive",
    outdir = results_path,
    meta_path = meta_path,
    sheet_name = sheet_name
    )

# 1.1.4 Filter (in) the Thy_vs_mut
uniqueRearr %>% 
  filter(!str_detect(`sample`, "Myc")) %>% 
  filter(!str_detect(`sample`, "WKM-CG1")) %>% 
  filter(!str_detect(`sample`, "WKM-CG1-P2mut")) %>% #View()
  plotBoxplotDotJitter(
    grouping = "Manuscript notation", 
    outname = "Thy_vs_mut", 
    ylab  = paste0("#clonotypes(",receptor_type,")"),
    subDir = "productive",
    outdir = results_path,
    meta_path = meta_path,
    sheet_name = sheet_name
  )

uniqueRearr %>% 
  filter(!str_detect(`sample`, "Myc")) %>% 
  filter(!str_detect(`sample`, "WKM-CG1")) %>% 
  filter(!str_detect(`sample`, "WKM-CG1-P2mut")) %>% #View()
  plotBoxplotDotJitterHz(
    grouping = "Manuscript notation", 
    outname = "Thy_vs_mut", 
    ylab  = paste0("#clonotypes(",receptor_type,")"),
    subDir = "productive",
    outdir = results_path,
    meta_path = meta_path,
    sheet_name = sheet_name
    )



