PATH=$/xtdisk/yangyg_group/zhoujy

paste reads.count $PATH/Smed/RNA-seq/DEG/Smed_length |awk '$1==$3'|cut -f1,2,4>i

awk 'NR==FNR{sum+=$2;next}{print $1,10^9*$2/(sum*$3)}' i i|head


for i in `cat 11`;do paste $i.count Smed_length |awk '$1==$3'|cut -f1,2,4 >$i.out; awk 'NR==FNR{sum+=$2;next}{print $1,10^9*$2/(sum*$3)}' $i.out $i.out >$i.rpkm;done

