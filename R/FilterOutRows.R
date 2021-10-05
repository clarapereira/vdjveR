filterOutRows <- function(df, string="TR",column = "allVHitsWithScore"){
  #'
  #' filter out rows from dataframe, df, which contains a string, string on specified column, column
  #' 
  #' @param df
  #' @param string
  #' @param column
  #' @example igh_germdj_mut <- igh_germdj %>% filterOutRows(string="TR", column = "allVHitsWithScore")
  #'
  df %>% filter(!str_detect(.data[[column]], "TR"))
}