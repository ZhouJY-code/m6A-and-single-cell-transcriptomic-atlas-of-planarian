library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(DoubletFinder)

siCTRL0d.data <- Read10X("*/siCTRL0d/cut-ctrl-0d/outs/filtered_feature_bc_matrix")
siCTRL6h.data <- Read10X("*/siCTRL6h/cut-ctrl-6h/outs/filtered_feature_bc_matrix")
siCTRL3d.data <- Read10X("*/siCTRL3d/cut-ctrl-3d/outs/filtered_feature_bc_matrix")
siCTRL7d.data <- Read10X("*/siCTRL7d/cut-ctrl-7d/outs/filtered_feature_bc_matrix")
siCTRL11d.data <- Read10X("*/siCTRL11d/cut-ctrl-11d/outs/filtered_feature_bc_matrix")

siWTAP0d.data <- Read10X("*/siWTAP0d/cut-kd-0d/outs/filtered_feature_bc_matrix")
siWTAP6h.data <- Read10X("*/siWTAP6h/cut-kd-6h/outs/filtered_feature_bc_matrix")
siWTAP3d.data <- Read10X("*/siWTAP3d/cut-kd-3d/outs/filtered_feature_bc_matrix")
siWTAP7d.data <- Read10X("*/siWTAP7d/cut-kd-7d/outs/filtered_feature_bc_matrix")
siWTAP11d.data <- Read10X("*/siWTAP11d/cut-kd-11d/outs/filtered_feature_bc_matrix")

colnames(x = siCTRL0d.data) <- paste('siCTRL0d', colnames(x = siCTRL0d.data), sep = '_')
colnames(x = siCTRL6h.data) <- paste('siCTRL6h', colnames(x = siCTRL6h.data), sep = '_')
colnames(x = siCTRL3d.data) <- paste('siCTRL3d', colnames(x = siCTRL3d.data), sep = '_')
colnames(x = siCTRL7d.data) <- paste('siCTRL7d', colnames(x = siCTRL7d.data), sep = '_')
colnames(x = siCTRL11d.data) <- paste('siCTRL11d', colnames(x = siCTRL11d.data), sep = '_')

colnames(x = siWTAP0d.data) <- paste('siWTAP0d', colnames(x = siWTAP0d.data), sep = '_')
colnames(x = siWTAP6h.data) <- paste('siWTAP6h', colnames(x = siWTAP6h.data), sep = '_')
colnames(x = siWTAP3d.data) <- paste('siWTAP3d', colnames(x = siWTAP3d.data), sep = '_')
colnames(x = siWTAP7d.data) <- paste('siWTAP7d', colnames(x = siWTAP7d.data), sep = '_')
colnames(x = siWTAP11d.data) <- paste('siWTAP11d', colnames(x = siWTAP11d.data), sep = '_')

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




siCTRL6h <- CreateSeuratObject(counts=siCTRL6h.data, project="siCTRL6h", min.cells = 15)
siCTRL6h[["percent.mt"]] <- PercentageFeatureSet(object = siCTRL6h, features = c("SMED30000702","SMED30001799","SMED30003686","SMED30010375","SMED30019959","SMED30031308"))
siCTRL6h <- subset(siCTRL6h, subset = nFeature_RNA > 500& nFeature_RNA < 6000 & percent.mt < 20)
siCTRL6h <- NormalizeData(siCTRL6h, normalization.method = "LogNormalize", scale.factor = 10000)
siCTRL6h <- FindVariableFeatures(siCTRL6h)
siCTRL6h@meta.data$tech <- "siCTRL6h"

combined.1 <- ScaleData(object = siCTRL6h, verbose = FALSE)
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
siCTRL6h_doublet  <- doubletFinder_v3(combined.6, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE)
save(siCTRL6h_doublet,file="siCTRL6h_doublet.Robj") 


siCTRL3d <- CreateSeuratObject(counts=siCTRL3d.data, project="siCTRL3d", min.cells = 15)
siCTRL3d[["percent.mt"]] <- PercentageFeatureSet(object = siCTRL3d, features = c("SMED30000702","SMED30001799","SMED30003686","SMED30010375","SMED30019959","SMED30031308"))
siCTRL3d <- subset(siCTRL3d, subset = nFeature_RNA > 500& nFeature_RNA < 6000 & percent.mt < 20)
siCTRL3d <- NormalizeData(siCTRL3d, normalization.method = "LogNormalize", scale.factor = 10000)
siCTRL3d <- FindVariableFeatures(siCTRL3d)
siCTRL3d@meta.data$tech <- "siCTRL3d"

combined.1 <- ScaleData(object = siCTRL3d, verbose = FALSE)
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
siCTRL3d_doublet  <- doubletFinder_v3(combined.6, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE)
save(siCTRL3d_doublet,file="siCTRL3d_doublet.Robj") 



siCTRL7d <- CreateSeuratObject(counts=siCTRL7d.data, project="siCTRL7d", min.cells = 15)
siCTRL7d[["percent.mt"]] <- PercentageFeatureSet(object = siCTRL7d, features = c("SMED30000702","SMED30001799","SMED30003686","SMED30010375","SMED30019959","SMED30031308"))
siCTRL7d <- subset(siCTRL7d, subset = nFeature_RNA > 500& nFeature_RNA < 6000 & percent.mt < 20)
siCTRL7d <- NormalizeData(siCTRL7d, normalization.method = "LogNormalize", scale.factor = 10000)
siCTRL7d <- FindVariableFeatures(siCTRL7d)
siCTRL7d@meta.data$tech <- "siCTRL7d"

combined.1 <- ScaleData(object = siCTRL7d, verbose = FALSE)
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
siCTRL7d_doublet  <- doubletFinder_v3(combined.6, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE)
save(siCTRL7d_doublet,file="siCTRL7d_doublet.Robj") 



siCTRL11d <- CreateSeuratObject(counts=siCTRL11d.data, project="siCTRL11d", min.cells = 15)
siCTRL11d[["percent.mt"]] <- PercentageFeatureSet(object = siCTRL11d, features = c("SMED30000702","SMED30001799","SMED30003686","SMED30010375","SMED30019959","SMED30031308"))
siCTRL11d <- subset(siCTRL11d, subset = nFeature_RNA > 500& nFeature_RNA < 6000 & percent.mt < 20)
siCTRL11d <- NormalizeData(siCTRL11d, normalization.method = "LogNormalize", scale.factor = 10000)
siCTRL11d <- FindVariableFeatures(siCTRL11d)
siCTRL11d@meta.data$tech <- "siCTRL11d"

combined.1 <- ScaleData(object = siCTRL11d, verbose = FALSE)
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
siCTRL11d_doublet  <- doubletFinder_v3(combined.6, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE)
save(siCTRL11d_doublet,file="siCTRL11d_doublet.Robj") 





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



siWTAP6h <- CreateSeuratObject(counts=siWTAP6h.data, project="siWTAP6h", min.cells = 15)
siWTAP6h[["percent.mt"]] <- PercentageFeatureSet(object = siWTAP6h, features = c("SMED30000702","SMED30001799","SMED30003686","SMED30010375","SMED30019959","SMED30031308"))
siWTAP6h <- subset(siWTAP6h, subset = nFeature_RNA > 500& nFeature_RNA < 6000 & percent.mt < 20)
siWTAP6h <- NormalizeData(siWTAP6h, normalization.method = "LogNormalize", scale.factor = 10000)
siWTAP6h <- FindVariableFeatures(siWTAP6h)
siWTAP6h@meta.data$tech <- "siWTAP6h"

combined.1 <- ScaleData(object = siWTAP6h, verbose = FALSE)
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
siWTAP6h_doublet  <- doubletFinder_v3(combined.6, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE)
save(siWTAP6h_doublet,file="siWTAP6h_doublet.Robj") 




siWTAP3d <- CreateSeuratObject(counts=siWTAP3d.data, project="siWTAP3d", min.cells = 15)
siWTAP3d[["percent.mt"]] <- PercentageFeatureSet(object = siWTAP3d, features = c("SMED30000702","SMED30001799","SMED30003686","SMED30010375","SMED30019959","SMED30031308"))
siWTAP3d <- subset(siWTAP3d, subset = nFeature_RNA > 500& nFeature_RNA < 6000 & percent.mt < 20)
siWTAP3d <- NormalizeData(siWTAP3d, normalization.method = "LogNormalize", scale.factor = 10000)
siWTAP3d <- FindVariableFeatures(siWTAP3d)
siWTAP3d@meta.data$tech <- "siWTAP3d"

combined.1 <- ScaleData(object = siWTAP3d, verbose = FALSE)
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
siWTAP3d_doublet  <- doubletFinder_v3(combined.6, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE)
save(siWTAP3d_doublet,file="siWTAP3d_doublet.Robj") 




siWTAP7d <- CreateSeuratObject(counts=siWTAP7d.data, project="siWTAP7d", min.cells = 15)
siWTAP7d[["percent.mt"]] <- PercentageFeatureSet(object = siWTAP7d, features = c("SMED30000702","SMED30001799","SMED30003686","SMED30010375","SMED30019959","SMED30031308"))
siWTAP7d <- subset(siWTAP7d, subset = nFeature_RNA > 500& nFeature_RNA < 6000 & percent.mt < 20)
siWTAP7d <- NormalizeData(siWTAP7d, normalization.method = "LogNormalize", scale.factor = 10000)
siWTAP7d <- FindVariableFeatures(siWTAP7d)
siWTAP7d@meta.data$tech <- "siWTAP7d"

combined.1 <- ScaleData(object = siWTAP7d, verbose = FALSE)
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
siWTAP7d_doublet  <- doubletFinder_v3(combined.6, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE)
save(siWTAP7d_doublet,file="siWTAP7d_doublet.Robj") 




siWTAP11d <- CreateSeuratObject(counts=siWTAP11d.data, project="siWTAP11d", min.cells = 15)
siWTAP11d[["percent.mt"]] <- PercentageFeatureSet(object = siWTAP11d, features = c("SMED30000702","SMED30001799","SMED30003686","SMED30010375","SMED30019959","SMED30031308"))
siWTAP11d <- subset(siWTAP11d, subset = nFeature_RNA > 500& nFeature_RNA < 6000 & percent.mt < 20)
siWTAP11d <- NormalizeData(siWTAP11d, normalization.method = "LogNormalize", scale.factor = 10000)
siWTAP11d <- FindVariableFeatures(siWTAP11d)
siWTAP11d@meta.data$tech <- "siWTAP11d"

combined.1 <- ScaleData(object = siWTAP11d, verbose = FALSE)
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
siWTAP11d_doublet  <- doubletFinder_v3(combined.6, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE)
save(siWTAP11d_doublet,file="siWTAP11d_doublet.Robj") 





merged<- merge(x=siCTRL0d, y=c(siCTRL6h,  siCTRL3d,  siCTRL7d,  siCTRL11d,  siWTAP0d,  siWTAP6h,  siWTAP3d,  siWTAP7d,  siWTAP11d) )

merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
save(merged,file="/xtdisk/yangyg_group/zhoujy/Smed/Seurat_1/merged_merge.Robj")    

merged <- RunPCA(merged, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)

save(merged,file="/xtdisk/yangyg_group/zhoujy/Smed/Seurat_1/merged_RunPCA.Robj")

pdf('ElbowPlot.pdf')
ElbowPlot(merged,ndims=100)
dev.off()

merged <- FindNeighbors(merged, reduction = "pca", dims = 1:20)
merged <- FindClusters(merged, resolution = 2, n.start = 10)
merged <- RunTSNE(merged, dims = 1:20,  nthreads = 4)
merged <- RunUMAP(merged, dims = 1:20,  nthreads = 4)

save(merged,file="/xtdisk/yangyg_group/zhoujy/Smed/Seurat_1/merged_findcluster.Robj")

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

write.table(all.markers,file="/xtdisk/yangyg_group/zhoujy/Smed/Seurat_1/allmarkers.txt")

