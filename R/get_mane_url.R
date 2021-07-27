#' @title get mane FTP urls
#' @description returns URLs from URLs that are valid (can be downloaded)
#' @export

get_mane_url <- function(){
  
  # helper functions
  is_url_valid <- function(url){
    HTTP_STATUS_OK <- 200
    hd <- httr::HEAD(url)
    status <- hd$all_headers[[1]]$status
    list(exists = status == HTTP_STATUS_OK, status = status)
  }
  
  # download lines
  ftp_dir <- 'https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/'
  lines <- readLines(ftp_dir)
  lines <- lines[grepl('release_',lines)]
  
  # check releases
  releases_raw <- unlist(lapply(strsplit(lines, split = '"'), function(x) x[2]))
  releases <- unlist(lapply(strsplit(gsub('\\/','',releases_raw), split = '_'), function(x) x[2]))
  last_release <- releases[length(releases)]
  
  # get file names
  site <- paste0('https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_',last_release,'/')
  RNA <- paste0('MANE.GRCh38.v',last_release,'.select_ensembl_rna.fna.gz')
  GFF <- paste0('MANE.GRCh38.v',last_release,'.select_ensembl_genomic.gff.gz')
  url_rna <- paste0(site, RNA)
  url_gff <- paste0(site, GFF)
  
  # check if file exists
  if (!is_url_valid(url_rna)$exists) stop('Something went wrong.. .rna url does not exist! ')
  if (!is_url_valid(url_gff)$exists) stop('Something went wrong.. .gff url does not exist! ')
  
  return(list(gff=url_gff, rna = url_rna))
  
}


