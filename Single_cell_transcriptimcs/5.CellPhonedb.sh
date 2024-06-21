
cellphonedb method statistical_analysis major_celltype.txt exp_matrix.txt --iterations=10 --threads=9
cellphonedb plot heatmap_plot --pvalues-path ./out/pvalues.txt --output-path ./ --count-name a.pdf major_celltype.txt

