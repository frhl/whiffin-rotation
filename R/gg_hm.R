#' @title create a colored heatmap
#' @description creates a colored barplot based on bonferroni correction.
#' 
#' @param df a data.frame with list_name and pvalue in data
#' @param bonf bonferroni correction threshold
#' @param colors colors used not significant, significant and bonferroni significant entries.
#' 
#' @family plotting
#' @export

gg_hm <- function(df, bonf, group, colors = c("firebrick2","white"), sig_labels = c('','*','**')){
  
  stopifnot(c('list_name','pvalue',group) %in% colnames(df))
  
  # log-pvalue
  df$logpvalue <- -log10(df$pvalue)
  
  # set significance thresholds
  df$significance <- sig_labels[1] #'Not significant'
  df$significance[df$pvalue < 0.05] <- sig_labels[2] #'Nominally significant'
  df$significance[df$pvalue < bonf] <- sig_labels[3] #'Bonferroni threshold'
  df$significance <- as.factor(df$significance)
  
  # generate plot
  plt <- ggplot(df, aes(y = reorder(list_name, logpvalue), x = df[[group]], fill = logpvalue, label = significance)) + 
    geom_tile(alpha = 0.9) + 
    geom_text() +
    scale_fill_gradient(low = colors[2], high = colors[1]) +
    xlab('Group') + ylab('GO Category') +
    theme_bw() 
  
  return(plt)
  
}