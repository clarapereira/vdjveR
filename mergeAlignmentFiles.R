#!/usr/bin/Rscript --vanilla

# Retrieves all files, merge into a single one and saves the data with all original columns and extra metadata columns 

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

samplename_sufix <- "_xcr"
all_productive_file_sufix <- ".ALL.productive.clones.txt"
all_file_sufix <- ".clonotypes.ALL.txt"

#outfilename <- "clonotypes_ALL_raw_wMeta_productive"


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

productive_clonotypes <- map(paste0(basedir,all_productive_file_sufix), data.table::fread) 
names(productive_clonotypes) <- samples_vector

all_clonotypes <- map(paste0(basedir,all_file_sufix), data.table::fread) 
names(all_clonotypes) <- samples_vector


nonproductive <-  map2( all_clonotypes, productive_clonotypes, getNonProductiveFromAll) 

# add class 
productive_clonotypes_mut <- map(productive_clonotypes, mutateWithString, string="productive")
nonproductive_mut <- map(nonproductive, mutateWithString, string="non_productive") 

# add sample names to rearrangements
productive_clonotypes_df <- productive_clonotypes_mut %>% 
  map2(samples_vector, mutateWithString, colname = sample_name_column_in_data) %>% 
  reduce(bind_rows) 
nonproductive_df <- nonproductive_mut %>% 
  map2(samples_vector, mutateWithString, colname = sample_name_column_in_data) %>% 
  reduce(bind_rows) 

# add metadata to rearrangements
productive_and_nonproductive_meta <- bind_rows(productive_clonotypes_df, nonproductive_df) %>% 
  left_join(
    metadata %>% rename(!! sample_name_column_in_data := sample_name_column_in_meta),
    by = sample_name_column_in_data 
  )


makeRatioColumns <- function(df, column_sufix = "all"){
  
  column_productive <- paste0(column_sufix, ".vdj+per_sample")
  column_nonproductive <- paste0(column_sufix, ".vdj-per_sample")
  `vdj+/vdj-` <- paste0(column_sufix, ".vdj+/vdj-")
  `vdj+/vdj+vdj-` <- paste0(column_sufix, ".vdj+/vdj+vdj-")
  # add productive/non-productive ratios:
  productive_per_sample_column <- df %>% 
    filter(rearrangement == "productive") %>% 
    group_by(sample) %>% 
    summarise(!! column_productive := n())
  nonproductive_per_sample_column <- df %>% 
    filter(rearrangement == "non_productive") %>% 
    group_by(sample) %>% 
    summarise(!! column_nonproductive := n())
  productive_ratio <- productive_per_sample_column %>% 
    left_join(
      nonproductive_per_sample_column,
      by = "sample"
    ) %>% 
    mutate(
      !! `vdj+/vdj-` := .data[[column_productive]]/.data[[column_nonproductive]],
      !! `vdj+/vdj+vdj-` := .data[[column_productive]]/(.data[[column_nonproductive]]+.data[[column_productive]])
    )
  
  return(productive_ratio)
}

all.productive_ratio <- productive_and_nonproductive_meta %>% makeRatioColumns()



# Export merged dataframe with all samples and metadata
allRatio.productive_and_nonproductive_meta <- productive_and_nonproductive_meta %>% 
  left_join(
    all.productive_ratio,
    by = "sample"
  ) 
allRatio.productive_and_nonproductive_meta %>% #nrow()
  write_delim(
    paste0(results_path, "/clonotypes_ALL_raw_wMeta_wRatios.all.tsv"),
    outfile, 
    delim = "\t"
  )

# filter for selected locus / receptor
TRB.productive_ratio <- productive_and_nonproductive_meta %>% 
  filter(str_detect(allVHitsWithScore, "TRB")) %>% 
  makeRatioColumns(column_sufix = "TRB") 
  
productive_and_nonproductive_meta.TRB <- allRatio.productive_and_nonproductive_meta %>% 
  filter(str_detect(allVHitsWithScore, "TRB")) %>% 
  left_join(
    TRB.productive_ratio,
    by = "sample"
  ) 

productive_and_nonproductive_meta.TRB %>% 
  write_delim(
    paste0(results_path, "/clonotypes_raw_wMeta_wRatios.TRB.tsv"),
    outfile, 
    delim = "\t"
  )

productive_meta.TRB <- productive_and_nonproductive_meta.TRB %>% 
  filter(rearrangement == "productive") 
productive_meta.TRB %>% 
  write_delim(
    paste0(results_path, "/clonotypes_raw_wMeta_wRatios.TRB.productive.tsv"),
    outfile, 
    delim = "\t"
  )
  
  
  
  
  
  
  
  
  
  