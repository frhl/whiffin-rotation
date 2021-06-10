#' @title get kozak strength
#' @description get kozak strength
#' @param x sequence
#' @export



get_utr_kozak_strength <- function(x){
  
  # THIS IS AN OLD FUNCTION. REDO!
  
  # strengths
  kozaks <- list(
    strong_1 = '(A|G)..ATGG',
    moderate_1 =   '(A|G)..ATG.',
    moderate_2 =   '...ATGG',
    weak_1 =   '...ATG.',
    weak_2 = 'ATG'
  )
  
  # get sequence around ATG
  orfs <- find_orfs(x)
  if (length(orfs) > 0){
    if (all(!is.na(orfs))){
      starts <- as.numeric(strsplit(names(orfs), split = '_')[[1]][1])
      check <- lapply(starts, function(i) ((i-5):(i+1)))
      check <- lapply(check, function(i) i[i>0])
      
      kozak_match <- lapply(check, function(is){
        
        seq_selected <- paste0(unlist(strsplit(x,split = ''))[is], collapse = '')
        # find kozak sequence
        lst <- lapply(kozaks, function(kozak) grepl(kozak, seq_selected))
        match <- min((1:5)[unlist(lst)])
        #names <- unlist(lapply(strsplit(names(kozaks), split = '_'), function(x) x[1]))
        names <- c(3,2,2,1,1)
        return(names[match])
        
      })
      
      names(kozak_match) <- starts
      kozak_match <- null_omit(kozak_match)
      return(kozak_match)
    }
  }
  return(NA)

}
