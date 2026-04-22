#####
#
# Astrangia poculata ICRT analysis and GO term
# Jake Warner
# April 22, 2026
#
#####

library(edgeR)

#read in the counts
counts <- read.table("AP_aligned.ap2.merged.txt", header=TRUE, sep='\t', row.names = 1)

colnames(counts) <- c(
  "DMSO_A",
  "DMSO_B",
  "DMSO_C",
  "ICRT_A",
  "ICRT_B",
  "SU_A",
  "SU_B"
)

#group them into replicates
group <- factor(c(1,1,1,2,2)) 
#setup the DGE object
y <- DGEList(counts=counts, group=group)
cpm.y <- cpm(y)

mean_lib <- mean(y$samples$lib.size)

y <- y[rowSums(cpm.y > ((5 / mean_lib)*1000000)) >=2,] 
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
z <- cpm(y, normalized.lib.size=TRUE)
write.table(z, file="normalized_cpm.txt", sep="\t", quote=F)

#DE testing is next:
y <- estimateCommonDisp(y)
sqrt(y$common.disp)
y <- estimateTagwiseDisp(y)
dim(y)
n_genes <-dim(y)[1]

de.DMSO.ICRT <- exactTest(y, pair=c(1,2))
tags.de.DMSO.ICRT <- topTags(de.DMSO.ICRT, n=(n_genes)) 
write.table(tags.de.DMSO.ICRT, file="tags.de.DMSO.ICRT", sep="\t", quote=F)

#get top tags
tags.de.DMSO.ICRT <- read.table(file = "tags.de.DMSO.ICRT", header=T, quote="", sep='\t',row.names = NULL)
colnames(tags.de.DMSO.ICRT)[1] <- "gene"

#get the union of the top tags
threshold_FC = 1.0
threshold_FDR = 0.05

top.tags.de.DMSO.ICRT <- tags.de.DMSO.ICRT[which(abs(tags.de.DMSO.ICRT$logFC)>= threshold_FC & tags.de.DMSO.ICRT$FDR<=threshold_FDR),]
write.table(top.tags.de.DMSO.ICRT, file="top.tags.de.DMSO.ICRT", sep="\t", quote=F)

##
## Volcano plot
##

tags.de.DMSO.ICRT <- read.table(file = "tags.de.DMSO.ICRT.Apoc2.1.txt", header=T, quote="", sep='\t')


annot<- read.table(file = "genomes/apoc2.1/apoculata_proteins_v2.1.renamed.new.totally_annotated.txt", sep='\t', quote="", header=T)
annot<- annot[c(1,9)]

skel_GOs <- "GO:0030198;|GO:0043062;|GO:0071711;|GO:0030199;|GO:0061448;|GO:0006029;|GO:0030166;|GO:0030203;|GO:0006024;|GO:0030204;|GO:0030206;|GO:0050654;|GO:0050651;|GO:0015015;|GO:0018146;|GO:0006816;|GO:0006811;|GO:0006812;|GO:0030001;|GO:0098660;|GO:0006874;|GO:0051480;|GO:0060429;|GO:0010631;|GO:0034109;|GO:0022407;|GO:0007155;|GO:0098609;|GO:0046903;|GO:1903530;|GO:0045055;|GO:1903530;"

annot$skel <- ifelse(
  apply(annot, 1, function(row) any(grepl(skel_GOs, row))),
  "yes",
  NA
)

threshold_FC = 1.5
threshold_FDR = 0.05

tags.de.DMSO.ICRT$diffexpressed <- "NO"
tags.de.DMSO.ICRT$diffexpressed[tags.de.DMSO.ICRT$logFC >= threshold_FC & tags.de.DMSO.ICRT$FDR<=threshold_FDR] <- "UP"
tags.de.DMSO.ICRT$diffexpressed[tags.de.DMSO.ICRT$logFC <= -(threshold_FC) & tags.de.DMSO.ICRT$FDR<=threshold_FDR] <- "DOWN"

tags.de.DMSO.ICRT.annotated <- merge(tags.de.DMSO.ICRT,annot, by.x='row.names', by.y='gene_ID', all.x=T)

#small<- tags.de.DMSO.ICRT.annotated[which(grepl("GO:0031214", tags.de.DMSO.ICRT.annotated$GO)),]

tags.de.DMSO.ICRT.annotated <- tags.de.DMSO.ICRT.annotated[rev(order(tags.de.DMSO.ICRT.annotated$skel)), ]
tags.de.DMSO.ICRT.annotated$skel <- factor(tags.de.DMSO.ICRT.annotated$skel, levels = c("yes"))

tags.de.DMSO.ICRT.annotated$label <- ifelse(tags.de.DMSO.ICRT.annotated$diffexpressed == "NO","not_significant", NA)
tags.de.DMSO.ICRT.annotated$label <- ifelse(tags.de.DMSO.ICRT.annotated$diffexpressed != "NO", "sig",tags.de.DMSO.ICRT.annotated$label)
tags.de.DMSO.ICRT.annotated$label <- ifelse((!is.na(tags.de.DMSO.ICRT.annotated$skel)) & tags.de.DMSO.ICRT.annotated$skel == "yes", "skel",tags.de.DMSO.ICRT.annotated$label)
tags.de.DMSO.ICRT.annotated$label <- ifelse(tags.de.DMSO.ICRT.annotated$diffexpressed == "NO","not_significant", tags.de.DMSO.ICRT.annotated$label)

tags.de.DMSO.ICRT.annotated$alpha <- ifelse(tags.de.DMSO.ICRT.annotated$label == "not_significant",0.5, NA)
tags.de.DMSO.ICRT.annotated$alpha <- ifelse(tags.de.DMSO.ICRT.annotated$label == "sig",1, tags.de.DMSO.ICRT.annotated$alpha)
tags.de.DMSO.ICRT.annotated$alpha <- ifelse(tags.de.DMSO.ICRT.annotated$label == "skel",1, tags.de.DMSO.ICRT.annotated$alpha)


p <- ggplot(data=tags.de.DMSO.ICRT.annotated, aes(x=logFC, y=-log10(FDR), col=label, alpha=alpha)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values=c("grey80","grey50","purple"))
p


##
## TopGO
##

#read in counts
meta <- read.table(file="tags.de.DMSO.ICRT.Apoc2.1.txt", sep = '\t',header=T)
annot<- read.table(file = "genomes/apoc2.1/apoculata_proteins_v2.1.renamed.new.totally_annotated.txt", sep='\t', quote="", header=T)

#this code chunk matches the GOs to the regen moduels and parses them
GO_Terms <- annot[c(1,9)]
colnames(GO_Terms) <- c("apoc_ID","GO")

library(plyr)
library(topGO)
meta_GO <- merge(meta, GO_Terms, by.x = "row.names",by.y = "apoc_ID", all.x = TRUE)

#get just the GO
meta_GO_Bkg <- meta_GO[c(1,6)]
#remove transcripts with no terms:
meta_GO_Bkg <- meta_GO_Bkg[complete.cases(meta_GO_Bkg),]
meta_GO_Bkg[meta_GO_Bkg==""]<-NA
meta_GO_Bkg[meta_GO_Bkg=="No_GO_Codes"]<-NA

#write it:
write.table(meta_GO_Bkg, file = "apoc_meta_GO_Bkg_parsed.txt", sep="\t", quote = F, row.names=F)

#this reads in and parses the GO terms
geneID2GO <- readMappings(file = "apoc_meta_GO_Bkg_parsed.txt", IDsep = "; ")

#this pulls the gene names from the GO_IDs, then makes a vector of genes matching black
geneNames <- names(geneID2GO)
head(geneNames)


myInterestingGenes <- meta[c(which(meta$logFC <= 1.5 & meta$FDR < 0.05)),]
myInterestingGenes <- row.names(myInterestingGenes)

#this creates a string with 0 to 1 for match to black
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)

#populate the list
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata

#run the Fisher test
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher

#output the results
allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 1000)
#cutoff the results at 0.05
allRes <- allRes[c(which(allRes$classicFisher < 0.05)),]
file = paste("GO_clust_ICRT_TOpGO_fisher_BP.txt",sep="")
write.table(allRes, file = file, sep='\t', quote=F, row.names = F)


library(forcats)
library(ggplot2)  
selected <- read.table(file="selected_skel_new_GO.txt", sep = '\t', header=T, quote = "")
selected$Pvalue <- -log(selected$classicFisher)
selected$Term <- factor(selected$Term, levels = selected$Term[order(selected$Pvalue)])

# Barplot
ggplot(selected, aes(x=Term, y=Significant, fill=Pvalue)) + 
  geom_bar(stat = "identity") +
  coord_flip() + 
  scale_fill_gradient(low = "blue", high = "#d73027", 
                      name="-Log10 p value" ) +
  theme_bw();
