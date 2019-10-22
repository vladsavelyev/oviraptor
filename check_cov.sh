BAM=$1
FAI=$2
if [ -z $FAI ] ; then
    FAI=/g/data/gx8/local/production/bcbio/genomes/Hsapiens/GRCh37/viral/gdc-viral.fa.fai
fi

PFX="${BAM%.*}"
OUT=$PFX-completeness.txt

# running mosdepth
mosdepth -t 32 $PFX $BAM -n --thresholds 1,5,25 --by <(awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $FAI) 

# output file header
echo -e "virus\tsize\treads\tdepth\t1x\t5x\t25x" > $OUT

# output file contents
paste \
    <(samtools idxstats $BAM | grep -v "\*" | sort -k1,1) \
    <(zcat $PFX.regions.bed.gz | sort -k1,1) \
    <(zgrep -v ^# $PFX.thresholds.bed.gz | sort -k1,1) \
    | awk 'BEGIN {FS="\t"} { print $1 FS $2 FS $3 FS $8 FS $14/$2 FS $15/$2 FS $16/$2}' \
    | sort -r -n -k3 \
    >> $OUT

# echo'ing output file name
echo $OUT
