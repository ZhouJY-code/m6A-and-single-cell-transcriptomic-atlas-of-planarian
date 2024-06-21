library(dplyr)
library(Matrix)
library(monocle)


load('data.Robj')
load('sub_cell.Robj')

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
pd <- new('AnnotatedDataFrame', data = Neoblast_neuronal@meta.data)
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
normal <- newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())

normal <- estimateSizeFactors(normal)
normal <- estimateDispersions(normal)
expressed_genes <- row.names(subset(fData(normal)))

save(normal,file="sample_normal.Robj")


###monoclemarker
clustering_DEG_genes <- differentialGeneTest(normal,fullModelFormulaStr = '~type',cores = 20)
save(clustering_DEG_genes,file="differentialGeneTest.Robj")

normal <- setOrderingFilter(normal, ordering_genes = clustering_DEG_genes)
#normal <- setOrderingFilter(normal, ordering_genes = head(clustering_DEG_genes[order(clustering_DEG_genes[,4]),],20))

###max_components
normal <- reduceDimension(normal, max_components = 2, method = 'DDRTree')
save(normal,file="reduceDimension.Robj")
##
normal <- orderCells(normal)
save(normal,file="orderCells.Robj")



pdf(file="trajectory.pdf")

plot_cell_trajectory(normal, color_by = "type")
plot_cell_trajectory(normal, color_by = "State")
plot_cell_trajectory(normal, color_by = "Pseudotime")
plot_cell_trajectory(normal, color_by = "tech")

dev.off()

my_pseudotime_de <- differentialGeneTest(normal,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 20)

save(my_pseudotime_de,file="my_pseudotime_de.Robj")

write.table(my_pseudotime_de,file="my_pseudotime_de.txt")

c<-subset(my_pseudotime_de, qval < 10^-20)

write.table(c[order(c[,4]),],file="my_pseudotime_de_qval.txt")

sig_gene_names <- row.names(subset(my_pseudotime_de, qval < 10^-20))

pdf(file="pseudotime_gene.pdf")
plot_pseudotime_heatmap(normal[sig_gene_names,], num_clusters = 3,cores = 10,use_gene_short_name = TRUE,show_rownames = TRUE,return_heatmap = TRUE)
dev.off()

gene<-cutree(cluster_gene$tree_row, k = 3)
write.table(gene,file="cutree.txt")

