library(Seurat)
library(Matrix)
library(dplyr)
library(methods)


#file<-read.table("pub_v6_smed_data.txt")
#mydata <- CreateSeuratObject(counts = file,min.cells = 0,min.features = 0,project = 'mydata_scRNAseq')
#meta<-read.csv("Planaria_Seurat_annot.csv",row.names=1)
#mydata<-AddMetaData(mydata, metadata=meta, col.name = NULL)
#sub_neural<-subset(mydata, subset=`final_Id`  %in% c("neural progenitors","ChAT neurons 1","ChAT neurons 2","GABA neurons","spp_11+ neurons","npp_18+ neurons","cav_1+ neurons","otf+ cells 1","otf+ cells 2"))
#allen_reference<-sub_neural
#allen_reference <- SCTransform(allen_reference,verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:50)
#save(allen_reference,file="allen_reference.Robj")

load('allen_reference.Robj')
load('merged_res6_rename.Robj')

sub_cell_neuronal<-subset(merged, idents="Neuronal")
sub_cell_neuronal <- SCTransform(sub_cell_neuronal,verbose = FALSE)

anchors <- FindTransferAnchors(reference = allen_reference, query = sub_cell_neuronal, normalization.method = "SCT")
save(anchors,file="anchors.Robj")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$final_Id, prediction.assay = TRUE, weight.reduction = sub_cell_neuronal[["pca"]])
sub_cell_neuronal[["predictions"]] <- predictions.assay
DefaultAssay(sub_cell_neuronal) <- "predictions"
save(sub_cell_neuronal,file="cell_type_pred.Robj")

data.table<- as.matrix(GetAssayData(object = sub_cell_neuronal, assay = "predictions"))

write.table(data.table,file="predictions.txt",sep="\t")

DefaultAssay(sub_cell_neuronal) <- "RNA"

pdf("ElbowPlot_PCA.pdf")
ElbowPlot(sub_cell_neuronal, ndims = 50, reduction = "pca")
dev.off()

sub_cell_neuronal <- RunPCA(sub_cell_neuronal, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
sub_cell_neuronal <- FindNeighbors(sub_cell_neuronal, reduction = "pca", dims = 1:10)
sub_cell_neuronal <- FindClusters(sub_cell_neuronal, resolution = 0.1, n.start = 15)
sub_cell_neuronal <- RunUMAP(sub_cell_neuronal, reduction = "pca", dims = 1:10)
sub_cell_neuronal <- RunTSNE(sub_cell_neuronal, reduction = "pca", dims = 1:10)


pdf('UMAP_PCA10_t.pdf')
DimPlot(sub_cell_neuronal, reduction="umap",pt.size = 0.5,label=T)
dev.off()

pdf('TSNE_PCA10_t.pdf')
DimPlot(sub_cell_neuronal, reduction="tsne",pt.size = 0.5,label=T)
dev.off()

DefaultAssay(sub_cell_neuronal) <- "predictions"
pdf('celltype_1.pdf')
FeaturePlot(sub_cell_neuronal, reduction="tsne",features = c("neural progenitors","ChAT neurons 1","ChAT neurons 2","GABA neurons"))
dev.off()
pdf('celltype_2.pdf')
FeaturePlot(sub_cell_neuronal,  reduction="tsne",features = c("spp-11+ neurons","npp-18+ neurons","cav-1+ neurons","otf+ cells 1"))
dev.off()
pdf('celltype_3.pdf')
FeaturePlot(sub_cell_neuronal,  reduction="tsne",features = "otf+ cells 2")
dev.off()

