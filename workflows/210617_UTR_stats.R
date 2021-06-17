# author: Frederik Lassen (21-06-09)
# description: 

devtools::load_all()
complexity <- fread('derived/tables/210615_MANE.v0.93.UTR_features.txt', sep = '\t')


# u5
x1 <- (data.frame(
  total_seq = nrow(complexity),
  u5_len_mean = mean(complexity$u5_len),
  u5_aug_total = sum(complexity$u5_AUG),
  u5_orf_total = sum(complexity$u5_ORF),
  u5_oorf_total = sum(complexity$u5_oORF_all | complexity$u5_oORF_altered_cds)
))

x2 <- (data.frame(
  total_seq = nrow(complexity),
  u3_len_mean = mean(complexity$u3_len),
  u3_aug_total = sum(complexity$u3_AUG),
  u3_orf_total = sum(complexity$u3_ORF)
))





