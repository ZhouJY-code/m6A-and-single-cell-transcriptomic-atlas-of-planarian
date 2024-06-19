## Volcano Plot ## deg数量需要修改
DEG<-read.table("E:/DESeq2_output.txt",header=T,row.names=1)
plot(DEG$log2FoldChange, -log10(DEG$pvalue), pch = 20, xlab = "log2(Fold Change)", ylab = "-log10(p value)", las = 1, main = "0h",xlim=c(-10,10),cex.lab=2,cex.axis=2,cex.main=5)
points(DEG[DEG[,2] >= 0.5849 & DEG[,5]<0.01,2], -log10(DEG[DEG[,2] >= 0.5849 & DEG[,5]<0.01,5]), pch=16, col="red");
points(DEG[DEG[,2] <=-0.5849 & DEG[,5]<0.01,2], -log10(DEG[DEG[,2] <=-0.5849 & DEG[,5]<0.01,5]), pch=16, col="green");


abline(v=0.5849,lty=8)
abline(v=-0.5849,lty=8)
abline(h=-log10(0.01),lty=8)
text(8,80, "n=624",cex=3)
text(-8,80, "n=261",cex=3)

## Volcano Plot ##
DEG<-read.table("E:/DESeq2_output.txt",header=T,row.names=1)
plot(DEG$log2FoldChange, -log10(DEG$pvalue), pch = 20, xlab = "log2(Fold Change)", ylab = "-log10(p value)", las = 1, main = "6h",xlim=c(-10,10),cex.lab=2,cex.axis=2,cex.main=5)
points(DEG[DEG[,2] >= 0.5849 & DEG[,5]<0.01,2], -log10(DEG[DEG[,2] >= 0.5849 & DEG[,5]<0.01,5]), pch=16, col="red");
points(DEG[DEG[,2] <=-0.5849 & DEG[,5]<0.01,2], -log10(DEG[DEG[,2] <=-0.5849 & DEG[,5]<0.01,5]), pch=16, col="green");


abline(v=0.5849,lty=8)
abline(v=-0.5849,lty=8)
abline(h=-log10(0.01),lty=8)
text(8,80, "n=363",cex=3)
text(-8,80, "n=290",cex=3)


## Volcano Plot ##
DEG<-read.table("E:/DESeq2_output.txt",header=T,row.names=1)
plot(DEG$log2FoldChange, -log10(DEG$pvalue), pch = 20, xlab = "log2(Fold Change)", ylab = "-log10(p value)", las = 1, main = "3d",xlim=c(-10,10),cex.lab=2,cex.axis=2,cex.main=5)
points(DEG[DEG[,2] >= 0.5849 & DEG[,5]<0.01,2], -log10(DEG[DEG[,2] >= 0.5849 & DEG[,5]<0.01,5]), pch=16, col="red");
points(DEG[DEG[,2] <=-0.5849 & DEG[,5]<0.01,2], -log10(DEG[DEG[,2] <=-0.5849 & DEG[,5]<0.01,5]), pch=16, col="green");


abline(v=0.5849,lty=8)
abline(v=-0.5849,lty=8)
abline(h=-log10(0.01),lty=8)
text(8,150, "n=1116",cex=3)
text(-8,150, "n=714",cex=3)


## Volcano Plot ##
DEG<-read.table("E:/DESeq2_output.txt",header=T,row.names=1)
plot(DEG$log2FoldChange, -log10(DEG$pvalue), pch = 20, xlab = "log2(Fold Change)", ylab = "-log10(p value)", las = 1, main = "7d",xlim=c(-10,10),cex.lab=2,cex.axis=2,cex.main=5)
points(DEG[DEG[,2] >= 0.5849 & DEG[,5]<0.01,2], -log10(DEG[DEG[,2] >= 0.5849 & DEG[,5]<0.01,5]), pch=16, col="red");
points(DEG[DEG[,2] <=-0.5849 & DEG[,5]<0.01,2], -log10(DEG[DEG[,2] <=-0.5849 & DEG[,5]<0.01,5]), pch=16, col="green");


abline(v=0.5849,lty=8)
abline(v=-0.5849,lty=8)
abline(h=-log10(0.01),lty=8)
text(8,150, "n=2374",cex=3)
text(-8,150, "n=3054",cex=3)



## Volcano Plot ##
DEG<-read.table("E:/DESeq2_output.txt",header=T,row.names=1)
plot(DEG$log2FoldChange, -log10(DEG$pvalue), pch = 20, xlab = "log2(Fold Change)", ylab = "-log10(p value)", las = 1, main = "11d",xlim=c(-10,10),cex.lab=2,cex.axis=2,cex.main=5)
points(DEG[DEG[,2] >= 0.5849 & DEG[,5]<0.01,2], -log10(DEG[DEG[,2] >= 0.5849 & DEG[,5]<0.01,5]), pch=16, col="red");
points(DEG[DEG[,2] <=-0.5849 & DEG[,5]<0.01,2], -log10(DEG[DEG[,2] <=-0.5849 & DEG[,5]<0.01,5]), pch=16, col="green");


abline(v=0.5849,lty=8)
abline(v=-0.5849,lty=8)
abline(h=-log10(0.01),lty=8)
text(8,150, "n=4050",cex=3)
text(-8,150, "n=5055",cex=3)
