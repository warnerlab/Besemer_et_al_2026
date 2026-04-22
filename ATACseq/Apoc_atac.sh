#####
#
# Astrangia poculata ATACseq alignment and quantification
# Jake Warner
# April 22, 2026
#
######



for fastq in Ap_12hpf_A_S10_L001_R1_001.fastq.gz ;
do
  cd ~/atac/
  #
  root=`basename $fastq _L001_R1_001.fastq.gz`;
  cat ${root}_L001_R1_001.fastq.gz ${root}_L002_R1_001.fastq.gz ${root}_L003_R1_001.fastq.gz ${root}_L004_R1_001.fastq.gz > ${root}_R1.fastq.gz
  cat ${root}_L001_R2_001.fastq.gz ${root}_L002_R2_001.fastq.gz ${root}_L003_R2_001.fastq.gz ${root}_L004_R2_001.fastq.gz > ${root}_R2.fastq.gz
  #
  fastqc ${root}_R1.fastq.gz &
  fastqc ${root}_R2.fastq.gz &
  #
  sh skewer_bwa.sh \
	${root}_R1.fastq.gz \
	${root}_R2.fastq.gz \
	${root} \
	/data/Apoculata/genomes/apoculata_v2.0.fasta \
	bwa_out \
	>> "${root}.log" 2>&1
  #
  cd bwa_out
  samtools sort -@ 38 -o ${root}.bwa.mem.sorted.bam ${root}.bwa.mem.bam
  samtools index ${root}.bwa.mem.sorted.bam

  alignmentSieve -b ${root}.bwa.mem.sorted.bam \
	--numberOfProcessors max \
	--minMappingQuality 5 \
	--maxFragmentLength 850 \
	--samFlagExclude 256 \
	--ATACshift \
	-o ${root}.sorted.filtered.bam \
	--filterMetrics ${root}.metrics.txt

  samtools sort -@ 38 -o ${root}.sorted.filtered.sorted.bam ${root}.sorted.filtered.bam
  rm ${root}.sorted.filtered.bam
  samtools index ${root}.sorted.filtered.sorted.bam

  bamCoverage --bam ${root}.sorted.filtered.sorted.bam  -o ${root}.sorted.filtered.coverage.bw \
	--numberOfProcessors max \
    --binSize 10 \
    --maxFragmentLength 850 \
    --extendReads 
done;

samtools view Ap_12hpf_A_S10.sorted.filtered.sorted.bam -L Apoc_mito.bed -o /dev/null -U Apoc_12hpfA.bam
samtools view Ap_12hpf_B-1_S11.sorted.filtered.sorted.bam -L Apoc_mito.bed -o /dev/null -U Apoc_12hpfB.bam
samtools view Ap_24Hpf_A_S12.sorted.filtered.sorted.bam -L Apoc_mito.bed -o /dev/null -U Apoc_24hpfA.bam
samtools view Ap_24Hpf_B-2_S13.sorted.filtered.sorted.bam -L Apoc_mito.bed -o /dev/null -U Apoc_24hpfB.bam
samtools view Ap_36Hpf_A_S14.sorted.filtered.sorted.bam -L Apoc_mito.bed -o /dev/null -U Apoc_36hpfA.bam
samtools view Ap_36Hpf_B_S15.sorted.filtered.sorted.bam  -L Apoc_mito.bed -o /dev/null -U Apoc_36hpfB.bam

samtools sort -@ 38 -n -o Apoc_12hpfA.name.sorted.bam Apoc_12hpfA.bam
samtools sort -@ 38 -n -o Apoc_12hpfB.name.sorted.bam Apoc_12hpfB.bam
Genrich -t Apoc_12hpfA.name.sorted.bam,Apoc_12hpfB.name.sorted.bam -o Ap_no_mito_12hpf.genrich.narrowPeak  -j -y -r &

samtools sort -@ 38 -n -o Apoc_24hpfA.name.sorted.bam Apoc_24hpfA.bam
samtools sort -@ 38 -n -o Apoc_24hpfB.name.sorted.bam Apoc_24hpfB.bam
Genrich  -t Apoc_24hpfA.name.sorted.bam,Apoc_24hpfB.name.sorted.bam -o Ap_no_mito_24hpf.genrich.narrowPeak  -j -y -r &

samtools sort -@ 38 -n -o Apoc_36hpfA.name.sorted.bam Apoc_36hpfA.bam
samtools sort -@ 38 -n -o Apoc_36hpfB.name.sorted.bam Apoc_36hpfB.bam
Genrich  -t Apoc_36hpfA.name.sorted.bam,Apoc_36hpfB.name.sorted.bam -o Ap_no_mito_36hpf.genrich.narrowPeak  -j -y -r &

##
## Annotate peaks
##

mkdir HOMER
cd HOMER

.//bin/annotatePeaks.pl \
../Ap_no_mito_12hpf.genrich.narrowPeak \
/mnt/data_stash/databases/omes/Apoculata/genomes/apoculata_v2.0.fasta \
-gff3 /mnt/data_stash/databases/omes/Apoculata/genomes/apoculata_v2.0.renamed.new.gff3 > Homer_12hpf_genrich_nomito_annotated_peaks.txt

wc -l Homer_12hpf_genrich_nomito_annotated_peaks.txt
6640
awk '{print $8}' Homer_12hpf_genrich_nomito_annotated_peaks.txt | grep -c intron
1196
awk '{print $8}' Homer_12hpf_genrich_nomito_annotated_peaks.txt | grep -c exon
645
awk '{print $8}' Homer_12hpf_genrich_nomito_annotated_peaks.txt | grep -c promoter
2809
awk '{print $8}' Homer_12hpf_genrich_nomito_annotated_peaks.txt | grep -c Intergenic
1650
awk '{print $8}' Homer_12hpf_genrich_nomito_annotated_peaks.txt | grep -c TTS
339

/home/warnerj/HOMER/.//bin/annotatePeaks.pl \
../Ap_no_mito_24hpf.genrich.narrowPeak \
/mnt/data_stash/databases/omes/Apoculata/genomes/apoculata_v2.0.fasta \
-gff3 /mnt/data_stash/databases/omes/Apoculata/genomes/apoculata_v2.0.renamed.new.gff3 > Homer_24hpf_genrich_nomito_annotated_peaks.txt

wc -l Homer_24hpf_genrich_nomito_annotated_peaks.txt
8702
awk '{print $8}' Homer_24hpf_genrich_nomito_annotated_peaks.txt | grep -c intron
1492
awk '{print $8}' Homer_24hpf_genrich_nomito_annotated_peaks.txt | grep -c exon
588
awk '{print $8}' Homer_24hpf_genrich_nomito_annotated_peaks.txt | grep -c promoter
4075
awk '{print $8}' Homer_24hpf_genrich_nomito_annotated_peaks.txt | grep -c Intergenic
2135
awk '{print $8}' Homer_24hpf_genrich_nomito_annotated_peaks.txt | grep -c TTS
411

/home/warnerj/HOMER/.//bin/annotatePeaks.pl \
../Ap_no_mito_36hpf.genrich.narrowPeak \
/mnt/data_stash/databases/omes/Apoculata/genomes/apoculata_v2.0.fasta \
-gff3 /mnt/data_stash/databases/omes/Apoculata/genomes/apoculata_v2.0.renamed.new.gff3 > Homer_36hpf_genrich_nomito_annotated_peaks.txt

wc -l Homer_36hpf_genrich_nomito_annotated_peaks.txt
9265
awk '{print $8}' Homer_36hpf_genrich_nomito_annotated_peaks.txt | grep -c intron
1562
awk '{print $8}' Homer_36hpf_genrich_nomito_annotated_peaks.txt | grep -c exon
603
awk '{print $8}' Homer_36hpf_genrich_nomito_annotated_peaks.txt | grep -c promoter
4336
awk '{print $8}' Homer_36hpf_genrich_nomito_annotated_peaks.txt | grep -c Intergenic
2328
awk '{print $8}' Homer_36hpf_genrich_nomito_annotated_peaks.txt | grep -c TTS
435

# Merge peaks

cat Ap_no_mito_12hpf.genrich.narrowPeak Ap_no_mito_24hpf.genrich.narrowPeak Ap_no_mito_36hpf.genrich.narrowPeak > Ap.all.genrich.nomito.narrowPeak
sort -k1,1 -k2,2n Ap.all.genrich.nomito.narrowPeak > Ap.all.genrich.nomito.sorted.narrowPeak
mergeBed -i Ap.all.genrich.nomito.sorted.narrowPeak > Ap.all.genrich.nomito.sorted.narrowPeak.merged

/home/warnerj/HOMER/.//bin/annotatePeaks.pl \
../Ap.all.genrich.nomito.sorted.narrowPeak.merged \
/mnt/data_stash/databases/omes/Apoculata/genomes/apoculata_v2.0.fasta \
-gff3 /mnt/data_stash/databases/omes/Apoculata/genomes/apoculata_v2.0.renamed.new.gff3 > Homer_Ap_ALL_genrich_noMito_annotated_peaks.txt

wc -l Homer_Ap_ALL_genrich_noMito_annotated_peaks.txt
11289
awk '{print $8}' Homer_Ap_ALL_genrich_noMito_annotated_peaks.txt | grep -c intron
2084
awk '{print $8}' Homer_Ap_ALL_genrich_noMito_annotated_peaks.txt | grep -c exon
930
awk '{print $8}' Homer_Ap_ALL_genrich_noMito_annotated_peaks.txt | grep -c promoter
4801
awk '{print $8}' Homer_Ap_ALL_genrich_noMito_annotated_peaks.txt | grep -c Intergenic
2877
awk '{print $8}' Homer_Ap_ALL_genrich_noMito_annotated_peaks.txt | grep -c TTS
596

##
## Transcription factor binding site annotation
## 

conda activate ananse
gimme scan Ap.all.genrich.nomito.sorted.narrowPeak.merged \
-g apoculata_v2.0 -p JASPAR2020_vertebrates -b > Apoc_no.mito.vert.motifs.bed
gimme motif2factors \
  --new-reference apoculata_v2.0 \
  --database JASPAR2020_vertebrates


