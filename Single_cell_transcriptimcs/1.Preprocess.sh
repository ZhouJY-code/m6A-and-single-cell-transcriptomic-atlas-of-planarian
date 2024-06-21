
export PATH=/software/biosoft/software/cellranger-3.0.0/:$PATH
data=cut-ctrl-0d

cd /home/Cellranger
cellranger count --id=cut-ctrl-0d --transcriptome=*/Smed_scRNA--fastqs=*/CellRanger/fastq_path --sample=cut-ctrl-0d --localmem=100 --localcores=30
