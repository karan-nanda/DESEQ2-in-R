library(ggplot2)
library(pheatmap)

normCount <- read.csv('Normal_Counts.csv', row.names = 1)
deSeqRes <- read.csv('Ordered_results.csv', row.names = 1)



deSeqRes$sig <- ifelse(deSeqRes$padj <= 0.05,'yes', 'no')


deSeqRes <- na.omit(deSeqRes)

ggplot(deSeqRes, aes(x = log10(baseMean), y=log2FoldChange, color = sig)) +
  geom_point()

