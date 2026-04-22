#######
#
# Astrangia poculata ATACseq timeseries analysis
# Jake Warner
# April 24, 2026
#
#######


library(edgeR)
counts <- read.table("Ap_timeseries.genrich.noMito.featurecounts.txt", header=TRUE, sep='\t', row.names = 1)
peaks_locations <- counts[c(1:5)]
counts <- counts[-c(1:5)]

colnames(counts) <- c(
  "12hpf_A",
  "12hpf_B",
  "24hpf_A",
  "24hpf_B",
  "36hpf_A",
  "36hpf_B")

group <- factor(c(1,1,2,2,3,3)) 
#setup the DGE object
y <- DGEList(counts=counts, group=group)
cpm.y <- cpm(y)


mean_lib <- mean(y$samples$lib.size)

#filter low OCRs
y <- y[rowSums(cpm.y > ((5 / mean_lib)*1000000)) >=1,] 
y <- calcNormFactors(y)
z <- cpm(y, normalized.lib.size=TRUE)

#setup for a PCA
#log2 transform
log_data = log(z+1,2)
transpose <- t(log_data)
transpose_df <- as.data.frame(transpose)
pca.data <- prcomp(transpose_df)
summary(pca.data) 

library(ggplot2)
scores = as.data.frame(pca.data$x) 

#plot it out
treatment <- factor(c("12hpf","12hpf","24hpf","24hpf","36hpf","36hpf"))
replicate <- factor(c(1,2,1,2,1,2))
p <- ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores), colour=factor(treatment), shape=factor(replicate))) + 
  geom_point(size=6) + 
  scale_fill_hue(l=40) +
  #scale_color_brewer(palette = "Purples") +
  coord_fixed(ratio=1, xlim=c(-80, 80), ylim=c(-80, 80)) + 
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  labs(title= "Astrangia development ATAC",color = "Timepoint", shape= "Replicate")
p

ggsave(p, filename = "Apoc_timeseries_PCA_all_norm.svg", height=5, width=6)


#DE testing is next:
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

dim(y)
n_peaks <-dim(y)[1]

#here's the de tests
de.12.24 <- exactTest(y, pair=c(1,2))
de.24.36 <- exactTest(y, pair=c(2,3))
de.12.36 <- exactTest(y, pair=c(1,3))
tags.de.12.24 <- topTags(de.12.24, n=(n_peaks))
tags.de.24.36 <- topTags(de.24.36, n=(n_peaks))
tags.de.12.36 <- topTags(de.12.36, n=(n_peaks))
write.table(tags.de.12.24, file="tags.peaks.de.12hpf.24hpf.allNorm.txt", sep="\t", quote=F)
write.table(tags.de.24.36, file="tags.peaks.de.24hpf.36hpf.allNorm.txt", sep="\t", quote=F)
write.table(tags.de.12.36, file="tags.peaks.de.12hpf.36hpf.allNorm.txt", sep="\t", quote=F)
tags.de.12.24 <- read.table(file="tags.peaks.de.12hpf.24hpf.allNorm.txt", sep="\t",row.names = 1)
tags.de.24.36 <- read.table(file="tags.peaks.de.24hpf.36hpf.allNorm.txt", sep="\t",row.names = 1)
tags.de.12.36 <- read.table(file="tags.peaks.de.12hpf.36hpf.allNorm.txt", sep="\t",row.names = 1)

#get coordinates
peaks <- read.table("Ap_timeseries.genrich.noMito.featurecounts.txt", header=TRUE, sep='\t', row.names = 1)

peaks <- peaks[c(1:3)]

de.12.24.peaks <- merge(peaks, tags.de.12.24, by.x='row.names', by.y='row.names', all.x=T)
row.names(de.12.24.peaks) <- de.12.24.peaks$Row.names
de.24.36.peaks <- merge(peaks, tags.de.24.36, by.x='row.names', by.y='row.names', all.x=T)
row.names(de.24.36.peaks) <- de.24.36.peaks$Row.names
de.12.36.peaks <- merge(peaks, tags.de.12.36, by.x='row.names', by.y='row.names', all.x=T)
row.names(de.12.36.peaks) <- de.12.36.peaks$Row.names

#testing output for chrom painting
de.12.24.peaks <- de.12.24.peaks[c(2:5,8)]
write.table(de.12.24.peaks, file="de.12.24.peaks.coords.allNorm.txt", sep="\t", quote=F, row.names =FALSE)
de.24.36.peaks <- de.24.36.peaks[c(2:5,8)]
write.table(de.24.36.peaks, file="de.24.36.peaks.coords.allNorm.txt", sep="\t", quote=F, row.names =FALSE)
de.12.36.peaks <- de.12.36.peaks[c(2:5,8)]
write.table(de.12.36.peaks, file="de.12.36.peaks.coords.allNorm.txt", sep="\t", quote=F, row.names =FALSE)


#need to remove chromosome_ from the coords
de.12.36.ideogram <- read.table("de.12.36.peaks.coords.allNorm.txt", sep = "\t", header = T, stringsAsFactors = F)
de.12.36.ideogram$Chr <- gsub("chromosome_", "", de.12.36.ideogram$Chr)
de.12.24.ideogram <- read.table("de.12.24.peaks.coords.allNorm.txt", sep = "\t", header = T, stringsAsFactors = F)
de.12.24.ideogram$Chr <- gsub("chromosome_", "", de.12.24.ideogram$Chr)
de.24.36.ideogram <- read.table("de.24.36.peaks.coords.allNorm.txt", sep = "\t", header = T, stringsAsFactors = F)
de.24.36.ideogram$Chr <- gsub("chromosome_", "", de.24.36.ideogram$Chr)

#how many signficantly DA?
volcano <- de.12.36.ideogram
volcano$diffexpressed[volcano$logFC > 0.6 & volcano$FDR < 0.05] <- "UP"
length(which(volcano$logFC > 0.6 & volcano$FDR < 0.05))
[1] 3189
volcano$diffexpressed[volcano$logFC < -0.6 & volcano$FDR < 0.05] <- "DOWN"
length(which(volcano$logFC < -0.6 & volcano$FDR < 0.05))
[1] 1159

#volcano Plot

myPalette <- c("#4575b4","#fc4d53") 
volcano_plot <- ggplot(data=volcano, aes(x=logFC, y=-log10(FDR), col=diffexpressed)) + 
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values = myPalette)
volcano_plot
ggsave(volcano_plot, filename = "Apoc_diff_access_volcano.12.36.allNorm.svg", height=6, width=5)

volcano <- de.12.24.ideogram
volcano$diffexpressed[volcano$logFC > 0.6 & volcano$FDR < 0.05] <- "UP"
length(which(volcano$logFC > 0.6 & volcano$FDR < 0.05))
[1] 2975
volcano$diffexpressed[volcano$logFC < -0.6 & volcano$FDR < 0.05] <- "DOWN"
length(which(volcano$logFC < -0.6 & volcano$FDR < 0.05))
[1] 1083

myPalette <- c("#4575b4","#fc4d53") 
volcano_plot <- ggplot(data=volcano, aes(x=logFC, y=-log10(FDR), col=diffexpressed)) + 
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values = myPalette)
volcano_plot
ggsave(volcano_plot, filename = "Apoc_diff_access_volcano.12.24.allNorm.svg", height=6, width=5)

volcano <- de.24.36.ideogram
volcano$diffexpressed[volcano$logFC > 0.6 & volcano$FDR < 0.05] <- "UP"
length(which(volcano$logFC > 0.6 & volcano$FDR < 0.05))
[1] 94
volcano$diffexpressed[volcano$logFC < -0.6 & volcano$FDR < 0.05] <- "DOWN"
length(which(volcano$logFC < -0.6 & volcano$FDR < 0.05))
[1] 37

myPalette <- c("#4575b4","#fc4d53") 
volcano_plot <- ggplot(data=volcano, aes(x=logFC, y=-log10(FDR), col=diffexpressed)) + 
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values = myPalette)
volcano_plot
ggsave(volcano_plot, filename = "Apoc_diff_access_volcano.24.36.allNorm.svg", height=6, width=5)

##
## chromosome painting
##

library(RIdeogram)
Apoc_karyotype <- read.table("../Apoc_chroms.txt", sep = "\t", header = T, stringsAsFactors = F)
ideogram(karyotype = Apoc_karyotype)
convertSVG("chromosome.svg", device = "png")
source("../ideaogram_jake.R")

## Ideogram

de.12.36.ideogram <- de.12.36.ideogram[c(1:4)]
colnames(de.12.36.ideogram) <- c("Chr", "Start", "End", "Value")

#enforce a scale limit
data <- sapply(de.12.36.ideogram$Value, function(x) ifelse(x>2, 2, x))
data <- sapply(data, function(x) ifelse(x<(-2), -2, x))
de.12.36.ideogram$Value <- data

ideogram_jake(karyotype = Apoc_karyotype, overlaid = de.12.36.ideogram, output = "chromosome_de.12.36.genrich.allNorm.svg")

##
## Making the OCR plots
##

library(ggplot2)

peak_annotations <- read.table(file="peak_annotations_tidy.txt", sep='\t', quote="",header=T)

to_plot<- peak_annotations[which(peak_annotations$location != "All"),]
to_plot<- to_plot[which(to_plot$condition == "hpf_12" | to_plot$condition ==  "hpf_24" | to_plot$condition == "hpf_36"),]

# Stacked
ggplot(to_plot, aes(fill=location, y=count, x=condition)) + 
  geom_bar(position="stack", stat="identity")

library(ggalluvial)
library(ggsci)
q <- ggplot(to_plot,
            aes(y = count, x = condition)) +
  geom_flow(aes(alluvium = location), alpha= .9, 
            lty = 2, fill = "white", color = "black",
            curve_type = "linear", 
            width = .5) +
  geom_col(aes(fill = location), width = .5, color = "black") +
  #scale_y_continuous(expand = c(0,0)) +
  scale_fill_brewer(palette = "Blues") +
  #scale_colour_brewer(palette = "Set1") +
  #scale_fill_manual(values=c("Exon"="violetred1","Intergenic"="seagreen1","Intron"= "yellow1","Promoter"="cyan")) +
  #pal_npg(palette = c("nrc"), alpha = 1)
  #scale_fill_npg() +
  cowplot::theme_minimal_hgrid()
q
ggsave(q, filename = "Apoc_timeseries_diff_access.svg", height=5, width=6)



