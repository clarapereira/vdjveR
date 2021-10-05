#!/usr/bin/Rscript --vanilla

# Retrieves all files, merge into a single one and saves the data with all original columns and extra metadata columns 
# performs equitability analysis, plots data.

# Note: 
# - the metadata file should contain a column with the names of samples, exactly as they are in the output file names of mixcr
# - each group of files for the same sample must be stored in a foldar with the same name as the sample
# - all those folders are inside the directory in the path "results_path"
# - meta_path is the path and name of the metadata file, as a .xlsx document
# - sheet_name is the name of the sheet containing the metadata


library(tidyverse)
source("./R/functions.R")


# inputs
meta_path <- "./zebrafish/meta/samples_zf_RNAseq_Mariana.xlsx"
sheet_name <- "main_editCP_missing_sample"
results_path <- "./zebrafish/results/"
sample_name_column_in_meta <- "sample_name"
sample_name_column_in_data <- "sample"

receptor_type <- "TRB"

samplename_sufix <- "_xcr"
filename_sufix <- ".ALL.productive.clones.txt"

outfilename <- "clonotypes_ALL_raw_wMeta_productive"


outfile <- paste0(results_path, "/", outfilename, ".tsv")


# Optional: add extra column with custum names to metadata
metadata <- readxl::read_xlsx(path = meta_path, sheet = sheet_name) %>% 
  distinct()

message(
  "Samples in metadata:\n",
  metadata %>% select(.data[[sample_name_column_in_meta]]) %>% distinct() 
)

# get all data
samples_vector <- metadata %>% filter(.data[[sample_name_column_in_meta]] != "NA") %>% select(.data[[sample_name_column_in_meta]]) %>% distinct() %>% na.omit() %>% pull()

basedir <- map(paste0(results_path, samples_vector, samplename_sufix, "/", samples_vector, samplename_sufix), paste0)
all_clonotypes <- map(paste0(basedir,filename_sufix), data.table::fread) 
names(all_clonotypes) <- samples_vector

# add metadata to rearrangements
all_clonotypes_meta <- all_clonotypes %>% 
  map2(samples_vector, mutateWithString, colname = sample_name_column_in_data) %>% 
  reduce(bind_rows) %>% 
  left_join(
    metadata %>% rename(!! sample_name_column_in_data := sample_name_column_in_meta),
    by = sample_name_column_in_data 
  )


# Export merged dataframe with all samples and metadata
all_clonotypes_meta %>% #nrow()
  write_delim(
    outfile, 
    delim = "\t"
  )


# =================================================================================================
# Up to here, the script is the same as mergeAlignmentFile.R
# =================================================================================================


# Make this modular (in independent script)




# =================================================================================================
# Add equitability for a particular receptor: 
# =================================================================================================

all_clonotypes_meta_w_equitability <- all_clonotypes_meta %>% 
  addEquitabilityPerSample(
    receptor_column = "allVHitsWithScore", 
    receptor_type = "TRB", 
    sample_column = "sample",
    outdir = results_path
  )

#all_clonotypes_meta_w_equitability %>% View()

# Equitability per sample per receptor type
receptorType_equitability <- all_clonotypes_meta_w_equitability  %>% 
  select(equitability_per_sample, sample) %>% 
  distinct()

# export TR dataframe w/ equitability
all_clonotypes_meta_w_equitability  %>% #View()
  write_delim(
    paste0(results_path, "/", outfilename, "_w_equitability.tsv"), 
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





