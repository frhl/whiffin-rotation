#' @title create a colored barplot
#' @description creates a colored barplot based on bonferroni correction.
#' 
#' @param df a data.frame with list_name and pvalue in data
#' @param bonf bonferroni correction threshold
#' @param colors colors used not significant, significant and bonferroni significant entries.
#' 
#' @family plotting
#' @export

gg_bar <- function(df, bonf, colors = c("firebrick2","tomato1","black")){
  
  stopifnot(c('list_name','pvalue') %in% colnames(df))
  
  # log-pvalue
  df$logpvalue <- -log10(df$pvalue)
  
  # set significance thresholds
  df$significance <- 'Not significant'
  df$significance[df$pvalue < 0.05] <- 'Nominally significant'
  df$significance[df$pvalue < bonf] <- 'Bonferroni threshold'
  df$significance <- as.factor(df$significance)
  
  # generate color palette
  names(colors) <- c('Bonferroni threshold','Nominally significant', 'Not significant')
  color_scale <- scale_fill_manual(name = "significance", values = colors)
  
  # generate plot
  plt <- ggplot(df, aes(y = reorder(list_name, logpvalue), x = logpvalue, fill = significance)) + 
    geom_bar(stat = 'identity', alpha = 0.9) +
    geom_vline(xintercept = -log10(0.05), color = "black", linetype = 'dashed') + 
    geom_vline(xintercept = -log10(bonf), color = "black", linetype = 'dashed') + 
    xlab('-log10(P-value)') + ylab('Experiment') +
    color_scale +
    theme_bw() 
  
  return(plt)
  
}