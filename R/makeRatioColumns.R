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