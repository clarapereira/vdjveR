#!/usr/bin/Rscript --vanilla

# Aim:
# Add equitability calculations, plot data.


# Usage: 
# Rscript doEquitabilityAnalysis.R


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
outfilename <- "clonotypes_raw_wMeta_wRatios.TRB.productive."

clonotypes <- data.table::fread(file)

# Optional: add extra column with custum names to metadata
metadata <- readxl::read_xlsx(path = meta_path, sheet = sheet_name) %>% 
  distinct()

# =================================================================================================
# Add equitability for a particular receptor: 
# =================================================================================================

clonotypes_w_equitability <- clonotypes %>% 
  addEquitabilityPerSample(
    receptor_column = "allVHitsWithScore", 
    receptor_type = "TRB", 
    sample_column = "sample",
    outdir = results_path
  )

#clonotypes_w_equitability %>% View()

# Equitability per sample per receptor type
receptorType_equitability <- clonotypes_w_equitability  %>% 
  select(equitability_per_sample, sample) %>% 
  distinct()

# export TR dataframe w/ equitability
clonotypes_w_equitability  %>% #View()
  write_delim(
    paste0(results_path, "/", outfilename, "w_equitability.tsv"), 
    delim = "\t"
  )



# =================================================================================================
# Plot the clonotypes and the equitability, independently:
# =================================================================================================
arrangement_type <- "productive"


# 1.1 plot number of unique clonotypes
# 1.1.2 Filter (in) the MYC samples: 
# 1.1.2.1 vertical 

receptorType_equitability %>% 
  filter(str_detect(`sample`, "Myc")) %>% 
  plotBoxplotDotJitter(
    y_max_lim = 1,
    grouping = "Manuscript notation", 
    ymeasurement = "equitability_per_sample",
    ylab  = paste0("Equitability(",receptor_type,")"),
    outname = "TransientLines", 
   # receptor_type = receptor_type#,
    #subDir = arrangement_type, 
    #outdir = results_path,
    #meta_path = meta_path,
    #sheet_name = sheet_name
  )
# 1.1.2.2 horizontal
receptorType_equitability %>% 
  filter(str_detect(`sample`, "Myc")) %>% 
  plotBoxplotDotJitterHz(
    y_max_lim = 1,
    grouping = "Manuscript notation", 
    ymeasurement = "equitability_per_sample",
    ylab  = paste0("Equitability(",receptor_type,")"),
    outname = "TransientLines",
    subDir = arrangement_type, 
    outdir = results_path,
    meta_path = meta_path,
    sheet_name = sheet_name
  )

# 1.1.3 Filter (in) the WKM_vs_stable
receptorType_equitability %>% 
  filter(!str_detect(`sample`, "Myc")) %>% 
  filter(!str_detect(`sample`, "Thy")) %>% 
  filter(!str_detect(`sample`, "WKM-CG1-P2mut")) %>% 
  plotBoxplotDotJitter(
    y_max_lim = 1,
    grouping = "Manuscript notation", 
    ymeasurement = "equitability_per_sample",
    outname = "WKM_vs_stable", 
    ylab  = paste0("Equitability(",receptor_type,")"),
    subDir = "productive",
    outdir = results_path,
    meta_path = meta_path,
    sheet_name = sheet_name
  )

receptorType_equitability %>% 
  filter(!str_detect(`sample`, "Myc")) %>% 
  filter(!str_detect(`sample`, "Thy")) %>% 
  filter(!str_detect(`sample`, "WKM-CG1-P2mut")) %>% 
  plotBoxplotDotJitterHz(
    y_max_lim = 1,
    grouping = "Manuscript notation", 
    ymeasurement = "equitability_per_sample",
    outname = "WKM_vs_stable", 
    ylab  = paste0("Equitability(",receptor_type,")"),
    subDir = "productive",
    outdir = results_path,
    meta_path = meta_path,
    sheet_name = sheet_name
    )

# 1.1.4 Filter (in) the Thy_vs_mut
receptorType_equitability %>% 
  filter(!str_detect(`sample`, "Myc")) %>% 
  filter(!str_detect(`sample`, "WKM-CG1")) %>% 
  filter(!str_detect(`sample`, "WKM-CG1-P2mut")) %>% #View()
  plotBoxplotDotJitter(
    y_max_lim = 1,
    grouping = "Manuscript notation", 
    ymeasurement = "equitability_per_sample",
    outname = "Thy_vs_mut", 
    ylab  = paste0("Equitability(",receptor_type,")"),
    subDir = "productive",
    outdir = results_path,
    meta_path = meta_path,
    sheet_name = sheet_name
  )

receptorType_equitability %>% 
  filter(!str_detect(`sample`, "Myc")) %>% 
  filter(!str_detect(`sample`, "WKM-CG1")) %>% 
  filter(!str_detect(`sample`, "WKM-CG1-P2mut")) %>% #View()
  plotBoxplotDotJitterHz(
    y_max_lim = 1,
    grouping = "Manuscript notation", 
    ymeasurement = "equitability_per_sample",
    outname = "Thy_vs_mut", 
    ylab  = paste0("Equitability(",receptor_type,")"),
    subDir = "productive",
    outdir = results_path,
    meta_path = meta_path,
    sheet_name = sheet_name
    )



# Add productive / non-productive ratio





