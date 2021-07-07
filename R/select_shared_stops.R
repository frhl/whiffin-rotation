



select_shared_stops <- function(x, mat, stop_shared){
  
  colnames(mat) <- c("Var1",'Var2')
  kozak <- stack(get_kozak_strength(x))
  colnames(kozak) <- c('str','Var1')
  kozak_mat <- merge(mat, kozak)
  kozak_mat$index <- 1:nrow(kozak_mat)
  
  keep <- unlist(lapply(unique(stop_shared), function(i){
    tmp_mat <- kozak_mat[kozak_mat$Var2 == i,]
    tmp_mat <- tmp_mat[tmp_mat$str == max(tmp_mat$str),]
    tmp_mat <- tmp_mat[tmp_mat$Var1 == min(tmp_mat$Var1),]
    return(tmp_mat$index)
  }))

  return(keep)
}