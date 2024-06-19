require(graphics)
require(grDevices)
library(gplots)
library(pheatmap)
a<-read.table("E:/expression1.txt",head=T,stringsAsFactors = F,row.names=1)
pheatmap(a,scale="row",cluster_rows = F,cluster_cols = F,show_colnames=T,show_rownames=T,legend=T,border_color="grey",color=colorpanel(60,"blue","white","red"),fontsize = 20,cellwidth = 60, cellheight = 60,gaps_row=c(1:2),cuttree_cols=4, gaps_col=c(1:4),cuttree_rows=2);
