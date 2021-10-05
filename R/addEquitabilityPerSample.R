addEquitabilityPerSample <- function(df, receptor_column = "allVHitsWithScore", receptor_type = "TRB", sample_column = "sample", outdir = results_path){
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
  
  # calculate equitability and add it to table
  df_w_equitability <-  df %>% 
    filter(str_detect(.data[[receptor_column]], receptor_type)) %>% 
    left_join(
      total_rearrangements_per_receptorType, 
      by = sample_column
    ) %>% 
    group_by(.data[[sample_column]]) %>% 
    mutate(
      total_cloneCount_per_sample = sum(cloneCount),
      clonotypeFraction_per_sample = cloneCount / total_cloneCount_per_sample,
      FlogF_per_sample = log(clonotypeFraction_per_sample)*clonotypeFraction_per_sample,
      sum_FlogF_per_sample = sum(FlogF_per_sample), # entropy
      # #total_unique_rearrangements = length(Groups),
      equitability_per_sample = -sum_FlogF_per_sample/log(total_unique_rearrangements_per_sample)
    ) %>% 
    ungroup()
  
  # Export equitability per sample per receptor type
  df_w_equitability  %>% 
    select(.data[[sample_column]], equitability_per_sample) %>% 
    mutate( receptor = receptor_type ) %>% 
    distinct() %>% 
    write_delim(paste0(outdir, "/equitability_persample_", receptor_type, ".tsv"), delim = "\t")
  
  return(df_w_equitability)
  
}