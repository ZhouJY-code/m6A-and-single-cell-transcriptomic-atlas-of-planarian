PATH=$/xtdisk/yangyg_group/zhoujy

cd $PATH/Smed/RNA-seq/macs2/0h/Ctrl-0h-1

macs2 callpeak -t ../../../MeRIP/7-Ctrl-0h-1-MeRIP_FKDL202579844-1a/map.bam  -c $PATH/Smed/RNA-seq/1.rawdata/7Ctrl0h-1_FRRB202356295-1a/map.bam -f BAMPE --nomodel  -g 6.3e8 -n m6a -q 0.001

cat m6a_peaks.xls|grep -v "^#"|sed '1,2d'|awk -v OFS="\t" '{print $1,$2,$3,"+"}' >ok.bed

cp $PATH/Smed/RNA-seq/macs2/distribution_plot_Smed.pl ./

cp $PATH/Smed/RNA-seq/macs2/Smed_ratio-distribution-mRNA ./

cp $PATH/Smed/RNA-seq/macs2/plot_tempfileR.txt ./

perl distribution_plot_Smed.pl

Rscript plot_tempfileR.txt

shuffleBed -i ok.bed  -g $PATH/Smed/RNA-seq/1.rawdata/index/Smed_length  > shuffle.bed
/pnas/yangyg_group/yangxin/software/Homer/bin/findMotifsGenome.pl ok.bed $PATH/Smed/RNA-seq/1.rawdata/index/smed_20140614.nt ./TOTAL/ -size given -p 1 -len 5,6,7 -rna -chopify -norevopp -cache 1000 -bg shuffle.bed

