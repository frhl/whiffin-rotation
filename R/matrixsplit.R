

matrixsplit <- function(mat, split, f = as.numeric, index = 3){
  
  return(apply(mat, 2, function(x) f(unlist(lapply(strsplit(x, split = split), function(x) x[index]))) ))
  
}