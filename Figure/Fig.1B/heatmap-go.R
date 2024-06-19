library(ggplot2)
input <-read.table("E:/go-heatmap.txt",header=T,sep="\t",stringsAsFactors=F)
rownames(input) <- input$Term
library(pheatmap)
pheatmap(as.matrix(-log10(input[,2:7])),color=colorpanel(60,"blue","white","red"),cellwidth =30, cellheight =30,show_rownames=T,border_color="gray",font_size=30,fontsize_row=20, fontsize_col=20,show_colnames=T,cluster_cols=F,cluster_rows=F,scale = "row",legend = T)



#pheatmap(as.matrix(-log10(input[,2:7])),color=colorpanel(100,"#008CC1","#FAF7F6","#CF0C31"),cellwidth =30, cellheight =30,show_rownames=T,font_size=30,fontsize_row=20, fontsize_col=20,show_colnames=T,cluster_cols=F,cluster_rows=F,scale = "row",legend = T,border=FALSE)