getNonProductiveFromAll <- function(df_all = igh_all, df_productive = igh_productive){
  #'
  #' @param df_all, the dataframe with all vdj alignments,
  #' @param df_productive, the dataframe with all vdj alignments excluded of out-of-frame and stop
  #' 
  #' @example igh_nonproductive <- igh_all[[1]] %>% getNonProductiveFromAll(igh_productive[[1]])
  #' 
  non_productive <- setdiff(
    df_all$aaSeqCDR3,
    df_productive$aaSeqCDR3
  )
  df_nonproductive <- df_all %>% #colnames()
    filter(  aaSeqCDR3 %in% non_productive )
  
  return(df_nonproductive)
}
