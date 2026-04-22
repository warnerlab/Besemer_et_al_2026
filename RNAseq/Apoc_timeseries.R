####
#
# Astrangia poculata developmental timeseries
# Jake Warner
# April 24, 2026
#
##### 

library(edgeR)

#
# apoculata timeseries
#

ap_counts <- read.table("Apoc_timeseries.featurecounts.txt", header=TRUE, sep='\t', row.names = 1)
ap_counts <- ap_counts[-c(1:5)]
group <- factor(c(1,1,2,2,3,3,4,4,5,5)) 
y <- DGEList(counts=ap_counts, group=group)
cpm.y <- cpm(y)
mean_lib <- mean(y$samples$lib.size)

y <- y[rowSums(cpm.y > ((5 / mean_lib)*1000000)) >=2,] 
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
z <- cpm(y, normalized.lib.size=TRUE)
#here's the cpm table for plotting etc.
write.table(z, file="apoc_normalized_cpm.txt", sep="\t", quote=F, row.names = T)

#pca
log_data = log(z+1,2)
transpose <- t(log_data)
transpose_df <- as.data.frame(transpose)
pca.data <- prcomp(transpose_df)
summary(pca.data) 

library(ggplot2)
scores = as.data.frame(pca.data$x) 
#plot it out
time <- factor(c("14hpf","14hpf","24hpf","24hpf","36hpf","36hpf","60hpf","60hpf","84hpf","84hpf"))
replicate <- factor(c(1,2,1,2,1,2,1,2,1,2))
p <- ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores), colour=factor(time), shape=factor(replicate))) + geom_point(size=6) + scale_fill_hue(l=40) + coord_fixed(ratio=1, xlim=c(-200, 150), ylim=c(-100, 100)) + 
  labs(title= "Astrangia poculata development",color = "Time Point", shape= "replicate")
p
ggsave(p, filename = "PCA_Apoc.pdf", height=8, width=12)


counts = read.table(file="apoc_normalized_cpm.txt", sep='\t', header=T)
master = log(counts+1,2)
master$avg_control_12hpf = (master$AP1A.APv2.0.sam +master$AP1B.APv2.0.sam)/2
master$avg_control_24hpf = (master$AP2A.APv2.0.sam +master$AP2B.APv2.0.sam)/2
master$avg_control_36hpf = (master$AP3A.APv2.0.sam +master$AP3B.APv2.0.sam)/2
master$avg_control_60hpf = (master$AP4A.APv2.0.sam +master$AP4B.APv2.0.sam)/2
master$avg_control_84hpf = (master$AP5A.APv2.0.sam +master$AP5B.APv2.0.sam)/2

master1 = merge(master,annotations, by.x='row.names',by.y='Apoc_ID', all.x=T)



write.table(master1,file="apoc_master_annotations.txt", sep='\t',row.names = F,quote=F)

rm(list = ls())
gc()


#
# Heatmap of GRN genes
#

rm(list = ls())
gc()

expression <- read.table(file="apoc_master_annotations.txt", sep='\t',quote="", header =T)
row.names(expression) <-  expression$Row.names
expression <-  expression[c(12:16)]

scaledata <- t(scale(t(expression))) # Centers and scales data.
scaledata <- scaledata[complete.cases(scaledata),]

grn <- read.table(file="Apoc_GRN_genes.txt", sep='\t',quote="", header =T)

scaledata_grn <- merge(scaledata,grn,by.x="row.names", by.y='Apoc_ID',all.y=T)

scaledata_grn <- scaledata_grn[complete.cases(scaledata_grn),]
row.names(scaledata_grn) <- scaledata_grn$gene_name
scaledata_grn <- scaledata_grn[c(2:6)]


library('viridis')
library(RColorBrewer)
colors = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(200)
library(WGCNA)
library("gplots")
mat<- as.matrix(scaledata_grn)
heatmap.2(x=mat)
pdf(file="plots/Apoc_GRN_heatmap.pdf", height = 40, width=6)
heatmap.2(x=mat, 
          Colv=FALSE, 
          dendrogram="none",
          scale="row",
          col=colors,
          trace="none",
          #sepcolor = "grey80",
          #colsep = 1:dim(scaledata_grn)[2],
          #sepwidth = c(0.01,0.01),
          density.info = "none",
          #main="my heatmap",
          ylab="Genes",
          xlab="Samples")
dev.off()
