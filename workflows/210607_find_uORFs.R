# author: Frederik Heymann Lassen, 07-June-2021
# description: find uORFs and various UTR complexity 
setwd('~/Projects/09_whiffin_rotation/whiffin-rotation/')

# load package
devtools::load_all()

# data that contains 5' UTRs, CDS and 3' UTRs
d <- fread('~/Projects/08_genesets/genesets/data/MANE/210705_MANE.GRCh38.v0.95.combined-table.txt', sep = '\t')
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
d2$u5_ORF <- unlist(lapply(d1$five_prime_UTR, function(x) length(get_orf(x, share_stops = F))))
d2$u5_oORF_all <-  unlist(lapply(1:nrow(d2), function(i) length(get_oorf(d1$five_prime_UTR[i], d1$CDS[i]))))
d2$u5_oORF_inframe <-  unlist(lapply(1:nrow(d2), function(i) length(get_oorf(d1$five_prime_UTR[i], d1$CDS[i], inframe = T))))
d2$u5_oORF_outframe <-  unlist(lapply(1:nrow(d2), function(i) length(get_oorf(d1$five_prime_UTR[i], d1$CDS[i], inframe = F))))
d2$u5_oORF_altered_cds <-  unlist(lapply(1:nrow(d2), function(i) length(get_altered_cds(d1$five_prime_UTR[i], d1$CDS[i]))))
d2$u5_oORF_kozak <-  unlist(lapply(1:nrow(d2), function(i) max(unlist(get_oorf_kozak(d1$five_prime_UTR[i], d1$CDS[i])))))
d2$u5_oORF_kozak[d2$u5_oORF_kozak == -Inf] <- 0 # the max on list operation returns warnings. fix here.
d2$u5_ORF_kozak <-  unlist(lapply(1:nrow(d2), function(i) max(unlist(get_orf_kozak(d1$five_prime_UTR[i])))))
d2$u5_ORF_kozak[d2$u5_ORF_kozak == -Inf] <- 0 # the max on list operation returns warnings. fix here.
d2$u5_ORF_cap_to_start <-  unlist(lapply(1:nrow(d2), function(i) get_cap_to_start_len(d1$five_prime_UTR[i])))
d2$u5_ORF_cds_proximity <-  unlist(lapply(1:nrow(d2), function(i) get_cds_proximity(d1$five_prime_UTR[i])))
d2$u5_ORF_leader_proximity <-  unlist(lapply(1:nrow(d2), function(i) get_leader_proximity(d1$five_prime_UTR[i])))
d2$u5_GC <- unlist(lapply(d1$five_prime_UTR, get_gc))

# 3' UTR
d2$u3_len <- unlist(lapply(d1$three_prime_UTR, function(x) nchar(x)))
d2$u3_AUG <- unlist(lapply(d1$three_prime_UTR, count_codon))
d2$u3_ORF <- unlist(lapply(d1$three_prime_UTR, function(x) length(get_orf(x))))
d2$u3_ORF_kozak <-  unlist(lapply(1:nrow(d2), function(i) max(unlist(get_orf_kozak(d1$three_prime_UTR[i])))))
d2$u3_ORF_kozak[d2$u3_ORF_kozak == -Inf] <- 0 # the max on list operation returns warnings. fix here.
d2$u3_GC <- unlist(lapply(d1$three_prime_UTR, get_gc))

warnings()
# write out table with UTR complexity features
d2 <- d2[!duplicated(d2),]
fwrite(d2, 'derived/tables/210707_MANE.v0.95.UTR_features.txt', sep = '\t')






