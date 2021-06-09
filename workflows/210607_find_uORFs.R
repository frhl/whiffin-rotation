
# process 5 and 3 prime UTRs to find complexity
source('R/find_orfs.R')
source('R/get_gc.R')
source('R/count_codon.R')
source('R/null_omit.R')
source('R/get_kozak_strength.R')
source('R/get_oorf.R')

library(data.table)

d <- fread('~/Projects/08_genesets/genesets/data/MANE/210607_MANE.GRCh38.v0.93.combined-table.txt', sep = '\t')
d$enstid <- unlist(lapply(strsplit(d$enstid_version, split = '\\.'), function(x) x[1]))

# get some UTRs stats
hist(unlist(lapply(d$seq[d$type == 'three_prime_UTR'], nchar)))
hist(unlist(lapply(d$seq[d$type == 'five_prime_UTR'], nchar)))

# get everything on one line
d1 <- do.call(rbind, lapply(unique(d$enstid_version), function(x){
  row <- d[d$enstid_version == x,]
  nrow <- row[,c(1,2,3,4,5)]
  nrow$five_prime_UTR <- row$seq[row$type == 'five_prime_UTR']
  nrow$CDS <- row$seq[row$type == 'CDS']
  nrow$three_prime_UTR <- row$seq[row$type == 'three_prime_UTR']
  return(nrow)
}))


# keep track of data here
d1 <- d1[!duplicated(d1),]
d2 <- d1
d2$five_prime_UTR <- NULL
d2$three_prime_UTR <- NULL
d2$CDS <- NULL

# 5' UTR
d2$u5_len <- unlist(lapply(d1$five_prime_UTR, function(x) nchar(x)))
d2$u5_AUG <- unlist(lapply(d1$five_prime_UTR, count_codon))
d2$u5_ORF <- unlist(lapply(d1$five_prime_UTR, function(x) length(find_orfs(x))))
d2$u5_oORF_all <-  unlist(lapply(1:nrow(d2), function(i) length(get_oorf(d1$five_prime_UTR[i], d1$CDS[i]))))
d2$u5_oORF_inframe <-  unlist(lapply(1:nrow(d2), function(i) length(get_oorf(d1$five_prime_UTR[i], d1$CDS[i], inframe = T))))
d2$u5_oORF_outframe <-  unlist(lapply(1:nrow(d2), function(i) length(get_oorf(d1$five_prime_UTR[i], d1$CDS[i], inframe = F))))
d2$u5_max_kozak <- unlist(lapply(d1$five_prime_UTR, function(x) max(unlist(get_utr_kozak_strength(x)))))
d2$u5_GC <- unlist(lapply(d1$five_prime_UTR, get_gc))

# 3' UTR
d2$u3_len <- unlist(lapply(d1$three_prime_UTR, function(x) nchar(x)))
d2$u3_AUG <- unlist(lapply(d1$three_prime_UTR, count_codon))
d2$u3_ORF <- unlist(lapply(d1$three_prime_UTR, function(x) length(find_orfs(x))))
d2$u3_max_kozak <- unlist(lapply(d1$three_prime_UTR, function(x) max(unlist(get_utr_kozak_strength(x)))))
d2$u3_GC <- unlist(lapply(d1$three_prime_UTR, get_gc))

d2 <- d2[!duplicated(d2),]
fwrite(d2, 'derived/tables/210610_MANE.v0.93.UTR_features.txt', sep = '\t')







# look at GC content
#plot(density(na.omit(d2$u5_GC)), main = 'GC content (%)', xlim = c(0.2,1), ylim = c(0, 4))
#lines(density(na.omit(gc_content[d$type == 'five_prime_UTR'])), col = 'red', legend = 'ga')
#lines(density(na.omit(gc_content[d$type == 'three_prime_UTR'])), col = 'blue')
#legend(0.91,4, legend=c("CDS", "5'-UTR", "3'-UTR"),
#       col=c("black", 'red',"blue"), lty=c(1,1,1), cex=1)

#d$five_prime_uorfs <- d$uorfs

#find_orfs(d[3,]$seq)
