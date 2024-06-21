require(graphics)
require(grDevices)
library(gplots)
library(pheatmap)
a<-read.table("E:/total-exp.txt",header=F,stringsAsFactors = F,row.names=1)

annotation_row = data.frame(
                    GeneClass = factor(rep(c("C1", "C2", "C3","C4"), c(4343,1410,941,640)))
                )
rownames(annotation_row) = rownames(a)
pheatmap(a,scale="row",cluster_rows = F,cluster_cols =F,show_colnames=T,show_rownames=F,legend=T,border_color="grey",color=colorpanel(60,"blue","white","red"),fontsize = 20,annotation_row = annotation_row);
