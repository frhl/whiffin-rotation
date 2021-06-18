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

hist(complexity$u5_ORF)

d1 <- as.data.frame(table(complexity$u5_ORF))
d1$Var1 <- as.numeric(as.character(d1$Var1))
d2a <- d1[d1$Var1 < 11,]
d2a$Var1 <- as.character(d2a$Var1)
d2b <- data.frame(Var1='10+', Freq  = sum(d1$Freq[d1$Var1 > 10]))
d2 <- rbind(d2a, d2b)



base <- complexity[complexity$u5_ORF != 0,]

max(base$u5_ORF)

d1 <- data.frame(x = 1:100, y = unlist(lapply(1:10, function(i) sum(complexity$u5_ORF == i))))#/nrow(base))))

ggplot(d1, aes(x = x, y = y)) +
  geom_point() +
  geom_line() + xlab('uORF number') + ylab('% of transcripts with at least one uORF') + ggtitle('n=6183')

#, type = 'l', xlab = 'uORF number', ylab = '% of transcripts with at least one uORF')

hist(complexity[complexity$u5_ORF != 0,]$u5_ORF, breaks = 30)



ggplot(d2, aes(x=Var1, y=Freq)) + geom_bar(stat='identity')




