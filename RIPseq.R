
Hfq <- read.table(file=(/Users/Nick/Desktop/Tobias/R/Hfq_0.5_min50_cov.txt), sep = '\t', header =F)
WT <- read.table(file=(/Users/Nick/Desktop/Tobias/R/WT0.5A.coverage.txt), sep = '\t', header =F)

names(Hfq) = c("A","B","C")
names(WT) = c("A","B","C")

Hfq$D <- WT$C[match(Hfq$B, WT$B)]

Hfq$E <- (Hfq$C/Hfq$D)

write.table(Hfq, file ="/Users/Nick/Desktop/test.txt", row.names = F, col.names =F)