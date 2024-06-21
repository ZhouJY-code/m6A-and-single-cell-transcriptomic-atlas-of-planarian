library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(DoubletFinder)

siCTRL0d.data <- Read10X("*/siCTRL0d/cut-ctrl-0d/outs/filtered_feature_bc_matrix")

siWTAP0d.data <- Read10X("*/siWTAP0d/cut-kd-0d/outs/filtered_feature_bc_matrix")

colnames(x = siCTRL0d.data) <- paste('siCTRL0d', colnames(x = siCTRL0d.data), sep = '_')

colnames(x = siWTAP0d.data) <- paste('siWTAP0d', colnames(x = siWTAP0d.data), sep = '_')

siCTRL0d <- CreateSeuratObject(counts=siCTRL0d.data, project="siCTRL0d", min.cells = 15)
siCTRL0d[["percent.mt"]] <- PercentageFeatureSet(object = siCTRL0d, features = c("SMED30000702","SMED30001799","SMED30003686","SMED30010375","SMED30019959","SMED30031308"))
siCTRL0d <- subset(siCTRL0d, subset = nFeature_RNA > 500& nFeature_RNA < 6000 & percent.mt < 20)
siCTRL0d <- NormalizeData(siCTRL0d, normalization.method = "LogNormalize", scale.factor = 10000)
siCTRL0d <- FindVariableFeatures(siCTRL0d)
siCTRL0d@meta.data$tech <- "siCTRL0d"

combined.1 <- ScaleData(object = siCTRL0d, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 20, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:20)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunUMAP(object = combined.4, reduction = "pca", dims = 1:20)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:20,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
mpK<-as.numeric(as.vector(combined.9$pK[which.max(combined.9$BCmetric)]))
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.01*length(combined.6@active.ident)) 
mpK<-as.numeric(as.vector(combined.9$pK[which.max(combined.9$BCmetric)]))
siCTRL0d_doublet  <- doubletFinder_v3(combined.6, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE)
save(siCTRL0d_doublet,file="siCTRL0d_doublet.Robj") 




siWTAP0d <- CreateSeuratObject(counts=siWTAP0d.data, project="siWTAP0d", min.cells = 15)
siWTAP0d[["percent.mt"]] <- PercentageFeatureSet(object = siWTAP0d, features = c("SMED30000702","SMED30001799","SMED30003686","SMED30010375","SMED30019959","SMED30031308"))
siWTAP0d <- subset(siWTAP0d, subset = nFeature_RNA > 500& nFeature_RNA < 6000 & percent.mt < 20)
siWTAP0d <- NormalizeData(siWTAP0d, normalization.method = "LogNormalize", scale.factor = 10000)
siWTAP0d <- FindVariableFeatures(siWTAP0d)
siWTAP0d@meta.data$tech <- "siWTAP0d"

combined.1 <- ScaleData(object = siWTAP0d, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 20, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:20)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunUMAP(object = combined.4, reduction = "pca", dims = 1:20)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:20,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0315*length(combined.6@active.ident)) 
mpK<-as.numeric(as.vector(combined.9$pK[which.max(combined.9$BCmetric)]))
siWTAP0d_doublet  <- doubletFinder_v3(combined.6, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE)
save(siWTAP0d_doublet,file="siWTAP0d_doublet.Robj") 



merged<- merge(x=siCTRL0d, y=c(siCTRL6h,  siCTRL3d,  siCTRL7d,  siCTRL11d,  siWTAP0d,  siWTAP6h,  siWTAP3d,  siWTAP7d,  siWTAP11d) )

merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
save(merged,file="merged_merge.Robj")    

merged <- RunPCA(merged, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)

save(merged,file="merged_RunPCA.Robj")

pdf('ElbowPlot.pdf')
ElbowPlot(merged,ndims=100)
dev.off()

merged <- FindNeighbors(merged, reduction = "pca", dims = 1:20)
merged <- FindClusters(merged, resolution = 2, n.start = 10)
merged <- RunTSNE(merged, dims = 1:20,  nthreads = 4)
merged <- RunUMAP(merged, dims = 1:20,  nthreads = 4)

save(merged,file="merged_findcluster.Robj")

pdf('TSNE.pdf')
DimPlot(merged, reduction = "tsne", pt.size = 0.1)
dev.off()

pdf('UMAP.pdf')
DimPlot(merged, reduction = "umap", pt.size = 0.1)
dev.off()

pdf('Group.pdf')
DimPlot(merged, reduction = "tsne", pt.size = 0.1, group.by="tech")
dev.off()

pdf('TSNE_lable.pdf')
DimPlot(merged, reduction = "tsne", pt.size = 0.1,label=T)
dev.off()

all.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.1)

write.table(all.markers,file="allmarkers.txt")

