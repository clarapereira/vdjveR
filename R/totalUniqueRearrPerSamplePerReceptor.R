totalUniqueRearrPerSamplePerReceptor <- function(df, receptor_column = "allVHitsWithScore", receptor_type = "TRB", sample_column = "sample", outdir){
  #' 
  #' Calculates equitability per sample and adds extra columns to the original table
  #' Takes a mixcr output dataframe, df, containing columns with sample identifiers (e.g., "sample"), and 
  #' columns containing the receptor type name (e.g., "allVHitsWithScore"), and, necessarily the column "cloneCount"
  #' 
  #'
  #' @param df
  #' @param receptor_column
  #' @param receptor_type
  #' @param sample_column
  #'
  #' @usage 
  
  
  # calculate total rearrangements in the full dataset and total rearr per sample: 
  total_rearrangements_per_receptorType <- df %>% #colnames()
    distinct() %>% 
    filter(str_detect(.data[[receptor_column]], receptor_type)) %>% 
    group_by(.data[[sample_column]]) %>% 
    summarise(total_unique_rearrangements_per_sample = n()) %>% 
    mutate(
      total_rearrangements = sum(total_unique_rearrangements_per_sample)
    )
  
  total_rearrangements_per_receptorType %>% 
  write_delim(paste0(outdir, "/unique_rearrangements_persample_", receptor_type, ".tsv"), delim = "\t")
  
  
  return(total_rearrangements_per_receptorType)
}