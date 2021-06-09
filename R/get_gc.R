

get_gc <- function(x){
 
 if (nchar(x) > 0){
   df <- as.data.frame(table(unlist(strsplit(x, split = ''))))
   colnames(df) <- c('base','freq') 
   return(sum(df$freq[df$base %in% c('G','C')])/sum(df$freq))
 } else {
   return(NA)
 }

 
}