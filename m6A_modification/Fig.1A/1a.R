require(graphics)
require(grDevices)
library(gplots)
library(pheatmap)
a<-read.table("E:/total.txt",header=F,stringsAsFactors = F,row.names=1)

annotation_row = data.frame(
                    GeneClass = factor(rep(c("C1", "C2", "C3","C4","C5","C6"), c(7677,3363,6853,6313,7334,3672)))
                )
rownames(annotation_row) = rownames(a)
pheatmap(a,scale="row",cluster_rows = F,cluster_cols =F,show_colnames=T,show_rownames=F,legend=T,border_color="grey",color=colorpanel(60,"blue","white","red"),fontsize = 20,annotation_row = annotation_row);


#pheatmap(a,scale="row",cluster_rows = F,cluster_cols =F,show_colnames=T,show_rownames=F,legend=T,border_color="grey",color=colorpanel(100,"#008CC1","#FAF7F6","#CF0C31"),fontsize = 20,annotation_row = annotation_row);