#!/bin/sh
#PBS -q core24
#PBS -l mem=20gb,walltime=500:00:00,nodes=1:ppn=5
#PBS -N index
#PBS -o index2.log
#PBS -e index2.err
#HSCHED -s empyro+perl+homo

cd /xtdisk/yangyg_group/sunbf/sunbf/Smed/RNA-seq/07123/Result-X101SC20052290-Z01-J003-B3-21/1.rawdata
for i in `cat 111`;do cd $i; bamCoverage -b sort.bam -o sort-c-1.bw --normalizeUsing CPM --binSize 1; cd ../; done
