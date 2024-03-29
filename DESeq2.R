library(DESeq2)
library(apeglm)

#loading the data

dat <- read.csv("airway_scaledcounts.csv", header = T, row.names  = 1)
info <- read.csv("airway_metadata.csv", header = T, row.names = 1)


deseqdf <- DESeqDataSetFromMatrix(dat, info, ~dex)

#Remove lowly expressed genes(row sum < 10)
keep <- rowSums(counts(deseqdf)) >= 10
deseqdf <- deseqdf[keep,]


#main dataframe

ddsDE <- DESeq((deseqdf))


#export normal read counts
normCounts <- counts(ddsDE, normalized = T)
write.csv(normCounts, 'Normal_Counts.csv')



#DESeq results
res <- results(ddsDE, alpha = 0.05) #Corrected p-value < 0.05

#output DESeq results
resOrdered <- res[order(res$padj),]
write.csv(resOrdered, 'Ordered_results.csv')

#summary(res)

#plotting logChange vs normalized counts
plotMA(ddsDE, ylim = c(-5,5))



