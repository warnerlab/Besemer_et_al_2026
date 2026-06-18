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


myInterestingGenes <- meta[c(which(meta$logFC <= -1.5 & meta$FDR < 0.05)),]
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
selected <- read.table(file="selected_secretory_new_GO.txt", sep = '\t', header=T, quote = "")
selected$Pvalue <- -log(selected$classicFisher)
selected$Term <- factor(selected$Term, levels = selected$Term[order(selected$Pvalue)])

  #lollipop plot
ggplot(selected, aes(x=Term, y=Significant, fill=Pvalue)) + 
  geom_bar(stat = "identity") +
  coord_flip() + 
  scale_fill_gradient(low = "blue", high = "#d73027", 
                      #limits = c( -4, 0),
                      #oob=squish, 
                      name="-Log10 p value" ) +
  theme_bw();
ggplot(selected, aes(x = reorder(Term, Significant), y = Significant)) +
  geom_segment(aes(xend = Term, y = 0, yend = Significant), 
               color = "grey50") +
  geom_point(aes(fill = Pvalue), 
             size = 4, 
             shape = 21, 
             color = "black", 
             stroke = 0.6) +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "#d73027",
                      name = "-Log10 p value") +
  labs(x = "GO Term", y = "Significant") +
  theme_bw()

