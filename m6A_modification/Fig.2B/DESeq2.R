export LD_LIBRARY_PATH=/software/biosoft/software/hdf51.10/lib:$LD_LIBRARY_PATH
export PATH=/software/biosoft/software/Rsoft/seurat3.1.2R/R3.6/bin:$PATH


#Ω¯»ÎR
library(DESeq2);

C1<-read.table("./Ctrl0h-1.count",head=F,stringsAsFactors=F);

C2<-read.table("./Ctrl0h-2.count",head=F,stringsAsFactors=F);

C3<-read.table("./Ctrl0h-3.count",head=F,stringsAsFactors=F);

T1<-read.table("./kd0h-1.count",head=F,stringsAsFactors=F);

T2<-read.table("./kd0h-2.count",head=F,stringsAsFactors=F);

T3<-read.table("./kd0h-3.count",head=F,stringsAsFactors=F);

data<-data.frame(C1[,2],C2[,2],C3[,2],T1[,2],T2[,2],T3[,2],stringsAsFactors=F);

rownames(data)<-c(C1[,1]);

colData<-data.frame(conditions=c(rep("C1",3),rep("T1",3)),row.names=colnames(data));

cds<-DESeqDataSetFromMatrix(data,design=~conditions,colData=colData);

cds<-DESeq(cds);

res <- results(cds)

resOrdered <- res[order(res$padj),]

plotMA(resOrdered);

addmargins( table( res_sig = resOrdered $padj < .1, res_sig = resOrdered $padj < .1 ) );

write.table (as.data.frame(resOrdered), file="DESeq2_output.txt",quote=F,sep="\t");

