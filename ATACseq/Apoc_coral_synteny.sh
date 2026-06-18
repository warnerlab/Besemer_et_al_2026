####
#
# Multi-way genome alignment and ATAC peak conservation
# June 18th 2026
# Jake Warner 
#
######

##
## Run cactus
##

vi evolverCorals.txt
#((((Mcapitata,Amillepora),(Spistillata,Gfascicularis)),Ofranksi),Apoculata);
Mcapitata /shared/ncbi_dataset/data/GCA_949126865.1_jaMonCapi2.1_genomic.chroms.fna 
Amillepora /shared/ncbi_dataset/data/GCF_013753865.1_Amil_v2.1_genomic.chroms.fna
Spistillata /shared/ncbi_dataset/data/GCA_964205215.1_jaStyPist1.1_genomic.chroms.fna
Gfascicularis /shared/ncbi_dataset/data/GCA_948470475.1_jaGalFasc40.1_genomic.chroms.fna 
Apoculata /shared/databases/omes/Apoculata/genomes/apoculata_v2.1.softmasked.fasta
Ofranksi Apoculata /shared/ncbi_dataset/dataGCA_964199315.1_jaOrbFran1.1_genomic.chroms.fna

vi run_coral_cactus.sh
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --partition=highmem
#SBATCH --mem=1800G
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=warnerj@uncw.edu
#SBATCH --mail-type=BEGIN,END,FAIL
singularity run /scripts_containers/cactus.sif \
cactus \
--maxMemory 1800G --maxCores 40 \
/projects/coral_synteny/jobstore \
/projects/coral_synteny/evolverCorals.txt \
/projects/coral_synteny/evolverCorals.hal

singularity run /scripts_containers/cactus.sif \
cactus-hal2maf --outType single --refGenome Apoculata /scratch/warnerj/coral_synteny/js evolverCorals.hal evolverCorals.Apoc_ref.maf

sbatch run_coral_cactus.sh

module load conda
conda activate /conda_envs/biopython

mkdir maf
mafSplit -byTarget -useFullSequenceName placeholder.bed maf/ evolverCorals.Apoc_ref.maf

##
## Phylofit and PhastCons
##

#get 4d sites from one chromosome
msa_view /projects/coral_synteny/maf/chromosome_1.maf --in-format MAF --4d \
 --features maf/chromosome_1.CDS.features.gtf \
 > 4d-codons.ss

#This will create a representation in the "sufficient statistics" (SS) format 
#of whole codons containing 4d sites. 4d sites (in the 3rd codon positions) can be extracted using msa_view.
msa_view 4d-codons.ss --in-format SS --out-format SS --tuple-size 1 > 4d-sites.ss

#phylofit to get non conserved
#tree came from Vaga et al. 2025 and pruned to target species
phyloFit --tree "((((Mcapitata,Amillepora),(Spistillata,Gfascicularis)),Ofranksi),Apoculata)" --msa-format SS --out-root nonconserved-4d 4d-sites.ss

#test phastCons
phastCons --target-coverage 0.25 --expected-length 12 \
 --rho 0.4 --estimate-rho mytrees --msa-format MAF \
 /projects/coral_synteny/maf/chromosome_1.maf nonconserved-4d.mod > scores.wig

# computed Rho =  (rho = 0.256564)

 phastCons --target-coverage 0.25 --expected-length 12 \
 --msa-format MAF /projects/coral_synteny/maf/chromosome_1.maf \
 --rho 0.256564 \
 --most-conserved most-cons.bed \
 mytrees.cons.mod,mytrees.noncons.mod > scores.wig

#loopity loop
for i in {1..14}; do \
  phastCons --target-coverage 0.25 --expected-length 12 \
  --msa-format MAF /projects/coral_synteny/maf/chromosome_$i.maf --rho 0.256564 \
  --most-conserved /projects/coral_synteny/maf/chromosome_$i.mostcons.bed \
  mytrees.cons.mod,mytrees.noncons.mod > /projects/coral_synteny/maf/chromosome_$i.phastcons.scores.bed.wig &
done

##
## Conservation of ATAC peaks
##

bedtools intersect -u \
-a Ap.2.1.all.genrich.nomito.sorted.narrowPeak.merged \
-b /projects/coral_synteny/Apoc_ref_mostcons.all.bed > \
 Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.bed

#all genes
wc -l Ap.2.1.all.genrich.nomito.sorted.narrowPeak.merged
wc -l Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.bed

annotatePeaks.pl \
Ap.2.1.all.genrich.nomito.sorted.narrowPeak.merged \
/databases/omes/Apoculata/genomes/apoculata_v2.1.fasta \
-gff3 /databases/omes/Apoculata/genomes/apoculata_v2.1.renamed.new.gff3 \
> Homer_Ap.all.genrich.sorted.narrowPeak.noMito.txt

wc -l Homer_Ap.all.genrich.sorted.narrowPeak.noMito.txt
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.txt | grep -c intron
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.txt | grep -c exon
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.txt | grep -c promoter
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.txt | grep -c Intergenic
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.txt | grep -c TTS


#most conserved genes
annotatePeaks.pl \
Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.bed \
/databases/omes/Apoculata/genomes/apoculata_v2.1.fasta \
-gff3 /databases/omes/Apoculata/genomes/apoculata_v2.1.renamed.new.gff3 \
> Homer_Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.txt

wc -l Homer_Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.txt
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.txt | grep -c intron
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.txt | grep -c exon
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.txt | grep -c promoter
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.txt | grep -c Intergenic
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.txt | grep -c TTS

#compare annotated to network genes
cat apoc_2.1_grn_ids.txt | while read line 
do
  grep --line-buffered "$line" Homer_Ap.all.genrich.sorted.narrowPeak.noMito.txt >> Homer_Ap.all.genrich.sorted.narrowPeak.noMito.GRN.txt
done
wc -l  Homer_Ap.all.genrich.sorted.narrowPeak.noMito.GRN.txt
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.GRN.txt | grep -c intron
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.GRN.txt | grep -c exon
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.GRN.txt | grep -c promoter
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.GRN.txt | grep -c Intergenic
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.GRN.txt | grep -c TTS

cat apoc_2.1_grn_ids.txt | while read line 
do
  grep --line-buffered "$line" Homer_Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.txt >> Homer_Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.GRN.txt
done
wc -l Homer_Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.GRN.txt
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.GRN.txt| grep -c intron
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.GRN.txt| grep -c exon
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.GRN.txt| grep -c promoter
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.GRN.txt| grep -c Intergenic
awk '{print $8}' Homer_Ap.all.genrich.sorted.narrowPeak.noMito.coral_syntenty.mostcons.GRN.txt| grep -c TTS
