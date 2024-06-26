export LD_LIBRARY_PATH=/software/biosoft/software/hdf51.10/lib:$LD_LIBRARY_PATH
export PATH=/software/biosoft/software/Rsoft/seurat3.1.2R/R3.6/bin:$PATH

library(topGO)
gid2GO <- readMappings(file ="./Smed_GOid1")
geneNames <-names(gid2GO)		    
geneTable <- read.table(file="1995-id",sep="\t",header=T,row.names=1)
GOI <- rownames(geneTable)                           
geneList<-factor(as.integer(geneNames %in% GOI))   
names(geneList)<-geneNames	 
GOenrich <- function(ontology="BP", allGenes, gene2GO)  
{
	GODat <- new("topGOdata", ontology=ontology, allGenes = allGenes,annot = annFUN.gene2GO, gene2GO=gene2GO)  
	resultCFis <- runTest(GODat, algorithm = "classic", statistic = "fisher") 
	allGO =usedGO(object = GODat)
	gtFis <- GenTable(GODat, classicFisher = resultCFis, orderBy= "classic", ranksOf = "classicFisher",  topNodes=length(allGO))
	GOgenes <- genesInTerm(GODat,gtFis$GO.ID) 
	sigGOI =sapply(names(GOgenes),function(x){paste(intersect(GOgenes[[x]],GOI),collapse=",")}) 
	sigGOI <- as.data.frame(sigGOI)      
	sigGOI <- cbind(sigGOI,"GO.ID"=rownames(sigGOI))   
	gtFis<-merge(gtFis,sigGOI,by.x="GO.ID",by.y="GO.ID")      
	gtFis[,6]<-as.numeric(gtFis[,6])   
	gtFis<-gtFis[order(gtFis[,6]),]   
	gtFis<-gtFis[gtFis[,6]<0.05,]     
	qv=p.adjust(gtFis[,6],method="fdr")
	FDR=formatC(qv,digits=4)
	temp<- gtFis[,7]
	gtFis[,7]<- FDR
	gtFis[,8]<- temp
	colnames(gtFis)<-c("GO.ID","Term","Annotated","Significant","Expected","classicFisher","FDR","Genelist")
	cbind("Class"=ontology,gtFis)	
}
header  <- c("Class","GO.ID","Term","Annotated","Significant","Expected","classicFisher","FDR","Genelist");
BPtable <- GOenrich(ontology="BP", allGenes=geneList, gene2GO=gid2GO);
colnames(BPtable)<-header  
write.table(BPtable,file="1995-id-topGO.xls",sep="\t",row.names=FALSE,col.names=TRUE)
