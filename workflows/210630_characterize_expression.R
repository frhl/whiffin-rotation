# 21-06-30: Get plots regarding expression across tissue

devtools::load_all()
library(cowplot)

# get protein / RNA 
expression <- fread('derived/tables/210609_prt_rna_numerical.txt', sep = '\t')

# plot protein  and RNA expression
p1 <- ggplot(expression, aes(x=rna, fill = tissue, color = tissue)) +
  geom_bar(show.legend = F) + 
  ylim(0,1000) +
  theme_bw() + 
  xlab('RNA abundance (log2 TPM)')

p2 <- ggplot(expression, aes(x=prt, fill = tissue, color = tissue)) +
  geom_bar() + 
  ylim(0,2500) +
  theme_bw() + 
  xlab('Protein abundance (log2 intensity)')

plot_grid(p1,p2,rel_widths = c(0.36,0.64))


# plot protein  and RNA expression
ggplot(expression, aes(x=prt, y = tissue, color = tissue)) +
  geom_boxplot() +
  theme_bw() + 
  xlab('RNA abundance (log2 TPM)')


plot_grid(p1,p2,rel_widths = c(0.36,0.64))


df$


library(rethinking)

