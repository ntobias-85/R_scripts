library("WGCNA")
library("DESeq2")
library("limma")
library("ape")

countdata <- read.table("Desktop/Tobias/R/Counts.txt", header=TRUE, row.names=1)
countdata <- as.matrix(countdata)

samples <- data.frame(row.names=c("WT_A","WT_B","hfq_KO_A","hfq_KO_B","hexA_KO_A","hexA_KO_B","hexA.hfq_KO_A","hexA.hfq_KO_B","hfq_comp_A","hfq_comp_B","WT_EV_hfq__A","WT_EV_hfq_B","PAL_KO_._EtOH_A","PAL_KO_._EtOH_B","PAL_KO_._CA_A","PAL_KO_._CA_B","plu4185_KO_A","plu4185_KO_B","plu4185_KO_comp_A","plu4185_KO_comp_B","WT_EV_plu4185_A","WT_EV_plu4185_B"),condition=as.factor(c(rep("WT",2),rep("hfq",2),rep("hexA",2),rep("hexAhfq",2),rep("hfq_comp",2),rep("WT_EV",2),rep("PAL_EtOH",2),rep("PAL_CA",2),rep("plu4185",2),rep("plu4185comp",2),rep("EV_plu4185",2))))

dds <- DESeqDataSetFromMatrix(countData=countdata, colData=samples, design=~condition)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
plot(sizeFactors(dds),colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))
logcounts <- log2( counts(dds, normalized=TRUE) + 1 )
plot(logcounts[,1], logcounts[,2], cex=.1)
rld <-rlog(dds)
plot(assay(rld)[,1], assay(rld)[,2], cex=.1)
normcounts <- counts(dds, normalized=TRUE)
write.table(normcounts, file="Desktop/Tobias/R/combined_noKO_norm_counts", sep="\t")

#Fix header

normdata <- read.table("Desktop/Tobias/R/normalized_counts.txt", header=TRUE, row.names=1)

#Once data is normalized
normdata <- as.vector(normdata)
#Takes the top 5000 genes - all in this case - and transposes
WGCNA_matrix = t(normdata[order(apply(normdata,1,mad), decreasing = T)[1:5000],])

#constructs a correlation matrix
s = abs(bicor(WGCNA_matrix))

# To pick the threshold and plot
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
      xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
      type='n', main = paste('Scale independence'));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
      labels=powers,cex=1,col='red'); abline(h=0.90,col='red')
#Here we select the threshold power as 'beta'
beta = 8
#create adjacency matrix
a = s^beta

exportNetworkToCytoscape(a,edgeFile="Desktop/Tobias/R/edge.txt", nodeFile="Desktop/Tobias/R/node.txt")