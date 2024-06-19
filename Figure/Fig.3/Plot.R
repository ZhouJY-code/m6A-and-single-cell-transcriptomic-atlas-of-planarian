
#3A
Idents(merged) <- factor(Idents(merged) , levels = c("Neoblast","Neuronal","Muscle","Parenchymal","Epidermal","Gut","Secretory"))
pdf('Umap_rename.pdf')
DimPlot(merged, reduction = "umap", pt.size = 0.1,label="F",cols=c("#CD378A","#816198","#B39775","#E4483D","#EE7079","#6BCCD9","#71B952"))
dev.off()
s
#3E
siCTRL<-SubsetData(merged, subset.name="orig.ident", accept.value= c("siCTRL0d","siCTRL6h","siCTRL3d","siCTRL7d","siCTRL11d"))
siWTAP<-SubsetData(merged, subset.name="orig.ident", accept.value= c("siWTAP0d","siWTAP6h","siWTAP3d","siWTAP7d","siWTAP11d"))
Idents(siCTRL) <- factor(Idents(siCTRL) , levels = c("Neoblast","Neuronal","Muscle","Parenchymal","Epidermal","Gut","Secretory"))
Idents(siWTAP) <- factor(Idents(siWTAP) , levels = c("Neoblast","Neuronal","Muscle","Parenchymal","Epidermal","Gut","Secretory"))

pdf('siWTAP_rename.pdf')
DimPlot(siWTAP, reduction = "umap", pt.size = 0.1,label="F",cols=c("#CD378A","#816198","#B39775","#E4483D","#EE7079","#6BCCD9","#71B952"))
dev.off()
pdf('siCTRL_rename.pdf')
DimPlot(siCTRL, reduction = "umap", pt.size = 0.1,label="F",cols=c("#CD378A","#816198","#B39775","#E4483D","#EE7079","#6BCCD9","#71B952"))
dev.off()

#3F
awk '$3>0&&$6<0.05{print $1}' Epidermal.txt > Epidermal_up.txt
awk '$3>0&&$6<0.05{print $1}' Gut.txt > Gut_up.txt
awk '$3>0&&$6<0.05{print $1}' Muscle.txt > Muscle_up.txt
awk '$3>0&&$6<0.05{print $1}' Neoblast.txt > Neoblast_up.txt
awk '$3>0&&$6<0.05{print $1}' Neuronal.txt > Neuronal_up.txt
awk '$3>0&&$6<0.05{print $1}' Parenchymal.txt > Parenchymal_up.txt
awk '$3>0&&$6<0.05{print $1}' Secretory.txt > Secretory_up.txt

awk '$3>0&&$6<0.05{print $1}' cut0d.txt > cut0d_up.txt
awk '$3>0&&$6<0.05{print $1}' cut6h.txt > cut6h_up.txt
awk '$3>0&&$6<0.05{print $1}' cut3d.txt > cut3d_up.txt
awk '$3>0&&$6<0.05{print $1}' cut7d.txt > cut7d_up.txt
awk '$3>0&&$6<0.05{print $1}' cut11d.txt > cut11d_up.txt

#3C
pdf('VlnPlot_neoblast.pdf',heigh=3)
VlnPlot(merged, features = "SMED30007406",pt.size = 0, combine = FALSE,cols=c("#CD378A","#816198","#B39775","#E4483D","#EE7079","#6BCCD9","#71B952"))
dev.off()

pdf('VlnPlot_neuronal.pdf',heigh=3)
VlnPlot(merged, features = "SMED30010096",pt.size = 0, combine = FALSE,cols=c("#CD378A","#816198","#B39775","#E4483D","#EE7079","#6BCCD9","#71B952"))
dev.off()

pdf('VlnPlot_muscle.pdf',heigh=3)
VlnPlot(merged, features = "SMED30003562",pt.size = 0, combine = FALSE,cols=c("#CD378A","#816198","#B39775","#E4483D","#EE7079","#6BCCD9","#71B952"))
dev.off()

pdf('VlnPlot_parenchymal.pdf',heigh=3)
VlnPlot(merged, features = "SMED30033522",pt.size = 0, combine = FALSE,cols=c("#CD378A","#816198","#B39775","#E4483D","#EE7079","#6BCCD9","#71B952"))
dev.off()

pdf('VlnPlot_epidermal.pdf',heigh=3)
VlnPlot(merged, features = "SMED30035904",pt.size = 0, combine = FALSE,cols=c("#CD378A","#816198","#B39775","#E4483D","#EE7079","#6BCCD9","#71B952"))
dev.off()

pdf('VlnPlot_gut.pdf',heigh=3)
VlnPlot(merged, features = "SMED30007290",pt.size = 0, combine = FALSE,cols=c("#CD378A","#816198","#B39775","#E4483D","#EE7079","#6BCCD9","#71B952"))
dev.off()

pdf('VlnPlot_secretory.pdf',heigh=3)
VlnPlot(merged, features = "SMED30001720",pt.size = 0, combine = FALSE,cols=c("#CD378A","#816198","#B39775","#E4483D","#EE7079","#6BCCD9","#71B952"))
dev.off()

