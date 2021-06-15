


gg_scatter_dens <- function(data, mapping, mapping.dens, title = '', subtitle = '', xlab = 'Normalized Log RNA expression', ylab = 'Normalized Log Protein Expression') {

  empty <- ggplot()+geom_point(aes(1,1), colour="white")+
    theme(axis.ticks=element_blank(), 
          panel.background=element_blank(), 
          axis.text.x=element_blank(), axis.text.y=element_blank(),           
          axis.title.x=element_blank(), axis.title.y=element_blank())
  
  plt <- ggplot(data, mapping) +
    geom_point(show.legend = F) +
    xlab(xlab) +
    ylab(ylab) +
    geom_hline(yintercept = 0, alpha = 0.2) +
    geom_vline(xintercept = 0, alpha = 0.2) +
    theme_bw()

  plt_rna <- ggplot(data, mapping.dens) +
    geom_density(alpha = 0.5, show.legend = T) +
    ggtitle(title, subtitle) +
    theme_bw() +
    xlab('') +
    theme(legend.position="top")
  
  plt_prt <- ggplot(data, mapping.dens) +
    geom_density(alpha = 0.5, show.legend = F) +
    coord_flip() +
    xlab('') +
    theme_bw()
  
  grid.arrange(plt_rna, empty, plt, plt_prt, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
  
  return(invisible(plt))
  
}

