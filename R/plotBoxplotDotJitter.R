plotBoxplotDotJitter <- function(df=total_rearrangements_receptorType , y_max_lim = 0 , grouping = "Groups", ymeasurement = "total_unique_rearrangements_per_sample", ylab  = paste0("#clonotypes(",receptor_type,")"), outname = "Groups", outdir = results_path, subDir = arrangement_type, meta_path, sheet_name){
  
  dir <- paste0(outdir,"/", subDir, "/")
  
  dir.create(file.path(dir), showWarnings = FALSE)
  
  #ylab  = paste0("#clonotypes(",receptor_type,")")
  
  
  
  d <- df %>% 
    left_join(
      metadata, 
      by =  c("sample"  = "sample_name")
    ) 
  ymaxlimit <- as.numeric(max(d %>% select(.data[[ymeasurement]])))
  
  set_ymax <- function(y_max_lim){
    if (y_max_lim == 1){
      1.05
    } 
    else 
      ymaxlimit+ymaxlimit*0.1
  }
  
  y_max_user <- set_ymax(y_max_lim)
  
  p <- d %>% #colnames()
    ggplot(aes(x=.data[[grouping]], y = .data[[ymeasurement]])) + 
    geom_boxplot(
      outlier.shape = NULL,
      outlier.size = NULL,
      outlier.stroke = 0
    ) +
    #geom_violin() +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, alpha = 0.5) +
    geom_jitter(aes(color = sample), shape=16, size = 3,  position=position_jitter(0.2), alpha = 0.8) +
    ylab(ylab) +
    ylim(0, y_max_user) +
    #theme(legend.position = "none") +
    theme_classic()  + 
    scale_color_manual(values = c26)# + # c25 is defined in the script colors.R
    #coord_flip()
  q <- p + 
    guides(color = FALSE) #+
    # theme(
    #   axis.text.x=element_text(angle = 45, hjust = 1)
    # )
  ggsave(
    paste0(dir,"plotBoxplotDotJitter", ylab ,"_",outname,".pdf"), 
    q,
    width = 2,
    height = 4)
  return(q)
}




