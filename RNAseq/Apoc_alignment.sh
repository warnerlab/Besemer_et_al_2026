##
#
# Alignment of Astrangia
#
##

hisat2-build /data/Apoculata/genomes/apoculata_v2.0.fasta

cd /mnt/data_stash/coral_embryos/data_dump/

#alignment

for fastq in AP*_1.fq.gz;
do
  root=`basename $fastq _1.fq.gz`;
  hisat2 -q -p 38 \
  --summary-file ${root}.APv2.hisat.aligned \
  -x /data/Apoculata/genomes/apoculata_v2.0.fasta \
  -1 ${root}_1.fq.gz \
  -2 ${root}_2.fq.gz \
  -S ${root}.APv2.0.sam
done
for sam in *.APv2.0.sam;
do
  root=`basename $sam .sam`;
  htseq-count --format=bam --stranded=no --type=mRNA --order=pos --idattr=ID \
  $sam \
  /data/Apoculata/genomes/apoculata_v2.0.renamed.gff3 > $root.counts.txt &
done
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' AP*counts.txt > AP_aligned.ap2.merged.txt

#feature counts

featureCounts -p -O -T 38 -t exon -g gene_id -a /data/Apoculata/genomes/apoculata_v2.0.renamed.gtf \
-o Apoc_timeseries.featurecounts.txt \
AP1A.APv2.0.sam \
AP1B.APv2.0.sam \
AP2A.APv2.0.sam \
AP2B.APv2.0.sam \
AP3A.APv2.0.sam \
AP3B.APv2.0.sam \
AP4A.APv2.0.sam \
AP4B.APv2.0.sam \
AP5A.APv2.0.sam \
AP5B.APv2.0.sam
