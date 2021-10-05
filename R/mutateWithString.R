mutateWithString <- function(df, colname = "rearrangement", string="non_productive"){
  #'
  #' add new column with a string
  #' @param df
  #' @param string
  #' @example igh_productive_mut <- igh_productive %>% mutateWithString(string="productive")
  #' 
  df_new <- df %>% 
    mutate( !!colname := string)
  return(df_new)
}