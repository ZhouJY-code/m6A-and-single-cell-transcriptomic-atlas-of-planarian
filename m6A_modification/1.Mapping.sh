PATH=/xtdisk/yangyg_group/zhoujy
smed="$PATH/Smed/RNA-seq/1.rawdata/index/smed_20140614.nt"
smed_length="$PATH/Smed/RNA-seq/1.rawdata/index/Smed_length"

cd $PATH/Smed/RNA-seq/1.rawdata

for i in `cat filename`;
do cd $i;


raw_1=${i}_1.fq.gz
raw_2=${i}_2.fq.gz


/software/biosoft/software/python/python2.7_2018_12/bin/cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o 1-trim_1.fq.gz -p 1-trim_2.fq.gz $raw_1 $raw_2;


java -Xmx4g -jar /pnas/yangyg_group/yangxin/software/new/Trimmomatic-0.33/trimmomatic-0.33.jar PE -phred33 ./1-trim_1.fq.gz ./1-trim_2.fq.gz ./2-trim_1.fq.gz ./unpaired_1.fq.gz ./2-trim_2.fq.gz ./unpaired_2.fq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:18;

rm 1-trim_1.fq.gz 1-trim_2.fq.gz;


bowtie2 -p 5 -N 1 -x $smed -1 2-trim_1.fq.gz -2 2-trim_2.fq.gz -S bowtie2.sam 2> log
samtools view -bS bowtie2.sam > map.bam

cat bowtie2.sam|grep -v '^@'|awk '$3!="*"'|cut -f1,3|awk '!a[$0]++'|cut -f2|sort|uniq -c|awk -vOFS="\t" '{print $2,$1}' > reads.count.tmp1

awk 'FNR==NR{a[$1]=$0;if(NR==1){NF--;gsub("[^ ]+","0");z=$0};next}{if($1 in a)print a[$1];else print $1"\t"z}' reads.count.tmp1 $smed_length > reads.count
###"-@ 30" should be deleted from samtools setence 
#samtools view -bS bowtie2.sam > bowtie2.bam
rm 2-trim_1.fq.gz
rm 2-trim_2.fq.gz

rm unpaired_1.fq.gz
rm unpaired_2.fq.gz

rm reads.count.tmp1
rm bowtie2.sam

cd ../;
done
