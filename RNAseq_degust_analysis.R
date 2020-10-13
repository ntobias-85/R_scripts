library(limma)
library(edgeR)

setwd(".")

count_cols <- c('Control_1_1','Control_1_2','Treatment_1_1','Treatment_1_2')
x<-read.delim('Count_data_merged.csv',skip=0, sep="	", check.names=FALSE, colClasses='character', na.strings=c())

x[,count_cols] <- apply(x[,count_cols], 2, function(v) as.numeric(v))     # Force numeric count columns
counts <- x[, count_cols]

row.names(x) = x[,1] # das war ich
row.names(counts) = x[,1]

keepMin <- apply(counts, 1, max) >= 0.0 # drop Genes. where no replicate features any mapped reads
keepCpm <- rowSums(cpm(counts)> 0.0) >= 0.0 # Keep only genes with cpm above x in at least y samples
keep <- keepMin & keepCpm
x <- x[keep,]
counts <- counts[keep,]
design <- matrix(c(c(1,1,0,0),c(0,0,1,1)), ncol=2, dimnames=list(c('Control_1_1','Control_1_2','Treatment_1_1','Treatment_1_2'),c('Control_1_','Treatment_1_')))


nf <- calcNormFactors(counts)
y<-voom(counts, design, plot=TRUE,lib.size=colSums(counts)*nf)

cont.matrix <- matrix(c(c(-1,1)), ncol=1, dimnames=list(c('Control_1_','Treatment_1_'),c('Treatment_1_')))

fit <- lmFit(y,design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

out <- topTable(fit2, n=Inf, sort.by='none')

out2 <- cbind(fit2$coef,
              out[, c('P.Value','adj.P.Val','AveExpr')],
              x[, c(c('Control_1_1','Control_1_2','Treatment_1_1','Treatment_1_2'))] )

write.csv(out2, file="DGE_results.csv", row.names=TRUE,na='',quote=FALSE) # changed row.names to TRUE, added quote=FALSE
