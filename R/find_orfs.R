





find_orfs <- function(x, start = 'ATG', stop = c('TAG', 'TAA', 'TGA')){
  
  # 
  seq <- unlist(strsplit(x, split = ''))
  starts <- strsplit()
  
  # get start codons
  split_start <- unlist(strsplit(x, split = start))
  split_start_n <- unlist(lapply(split_start, nchar))+3 
  starts <- cumsum(split_start_n) # start codon = i-2 : i .. e.g. seq[(1610-2):1610]
  
  # get stop codons
  split_stop <- unlist(strsplit(x, split = stop))
  split_stop_n <- unlist(lapply(split_stop, nchar))+3 
  stops <- cumsum(split_stop_n) # start codon = i-2 : i .. e.g. seq[(1927-2):1927]
  
  # find matching stop-start codons (ORFs)
  for (cur_start in starts){
    
    # look through stops 
    valid_stops <- stops[stops > cur_start]
    inframe_stops_bool <- (valid_stops-cur_start) %% 3 == 0
    inframe_stops <- valid_stops[inframe_stops_bool]
    
    # assertions
    stopifnot(paste0(seq[(cur_start-2):(cur_start)], collapse = '') == 'ATG')
    
    
    
    
    
    #for (cur_stop in inframe_stops){
    #  
    #  seq[(cur_start-2):(cur_stop)]
    #  
    #  stopifnot(length(seq[(cur_start-2):(cur_stop)]) %% 3 == 0)
    #  
    #}
    
    
    
    
    
  }
  
}
  
  
  