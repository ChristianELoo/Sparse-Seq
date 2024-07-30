cd /mnt/isilon/zhou_lab/projects/20211119_5hmC_Project
for in_bam in tmp/bam_bismarkbt2/GSM320718{1..6}.bam; do
  sname=$(basename ${in_bam} .bam)
  samtools idxstats "${in_bam}" | awk '{n+=$3} END {print n;}' >tmp/$sname.readcnt
done

mkdir -p tmp/subBam

for seed in {2..100}; do
  cat <<EOF | pbsgen -days 20 -name downsample_rep${seed} -submit
for in_bam in tmp/bam_bismarkbt2/GSM320718{1..6}.bam; do
  sname=\$(basename \${in_bam} .bam)
  for power in {6..18}; do 
    n_reads=\$(((2**\$power)/3*2)); ## here the downsampled read number is scaled to match read length
    echo 2**\$power \$n_reads
    /mnt/isilon/zhoulab/labtools/Rutils/wzsampleNumber.R \$(cat tmp/\$sname.readcnt) \$n_reads | awk 'NR==FNR{a[\$1]}NR!=FNR{if(\$0~/^@/) {print;} else {k+=1; if(k in a) {print;}}}' - <(samtools view -h \$in_bam) | samtools view -b >tmp/subBam/\${sname}_down_2pow\${power}_${seed}.bam
    /home/zhouw3/zhoulab/labsoftware/bismark/default/bismark_methylation_extractor --no_overlap --gzip tmp/subBam/\${sname}_down_2pow\${power}_${seed}.bam -o tmp/subBam/
    rm -f tmp/subBam/\${sname}_down_2pow\${power}_${seed}.bam
  done
done
EOF
done

# header: filename, power, seed, gsm, nreads,mCG,uCG,mCHG,uCHG,mCHH,uCHH;
for f in tmp/subBam/*splitting_report.txt; do awk '/Total number of methylation call strings processed:/{match($1, "Total number of methylation call strings processed: ([0-9]*)", aa); nreads=aa[1];}/Total methylated.*CpG context:/{mCG=$2;}/Total C to T conversions in CpG context/{uCG=$2;}/Total methylated.*CHG context:/{mCHG=$2;}/Total C to T conversions in CHG/{uCHG=$2;}/Total methylated.*CHH context:/{mCHH=$2;}/Total C to T conversions in CHH/{uCHH=$2;}END{match(FILENAME, "([^/]*)_down_2pow([0-9]*)_([0-9]*)", a); print FILENAME,a[2],a[3],a[1],nreads,mCG,uCG,mCHG,uCHG,mCHH,uCHH;}' $f; done >stats/subBam_counts_20230219.tsv