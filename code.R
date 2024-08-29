### Two datasets were used in this study. Count and matrix files can be downloaded directly from the NCBI GEO using the reference for each dataset.

###########################################################################################
##### Dataset 1: GSE57148
###########################################################################################

# load packages and sources file with function definitions.
library(tidyverse)
source(functions.R)

dataset <- read.csv("GSE57148/GSE57148_COPD_FPKM_Normalized.txt", sep = "\t")
rownames(dataset) <- dataset$GeneName
dataset$GeneName <- NULL

datasetmeta <- read.csv("GSE57148/SraRunTable.txt")
datasetmeta$sampleID <- colnames(dataset)

tpm <- fpkmToTpm(dataset)



# Dataset was only deposited as FPKM in GEO, therefore we used limma to do differential expression on the dataset.

datasetl <- log(dataset+0.1, base = 2)

library(limma)

### generate design matrix
predesign <- datasetmeta[, c("sampleID","disease_state")] %>% column_to_rownames("sampleID")
predesign$Normal <- 0
predesign$COPD <- 0
predesign[predesign$disease_state=="Normal",]$Normal <- 1
predesign[predesign$disease_state=="COPD",]$COPD <- 1
design <- predesign[, c("Normal", "COPD")]

#fit model
fit <- lmFit(datasetl, design = design)
cont.matrix <- makeContrasts(COPDvsNormal=COPD-Normal, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, trend = TRUE)

topTable(fit2, adjust.method = "BH")

results <- topTable(fit2, adjust.method = "BH", number = nrow(dataset))

## GENERATE VOLCANO Plot used in figure 3:
# load gags enzymes list
gags <- read.csv("data/gags.csv", header = T, stringsAsFactors = F, sep = ";")

results$color <- "other"
gags1 <- results[rownames(results)%in%gags$Gene,]

for (i in unique(gags$GAG)){
results[rownames(results)%in%gags[gags$GAG==i,]$Gene,]$color <- i
}
unique(results$color)
results[,1:6] <- as.data.frame(apply(results[,1:6], 2, as.numeric))
res <- results[results$color=="other",]
gags1 <- results[results$color!="other",]
gags1$label <- ""
gags1[abs(gags1$logFC)>=0.3 & -log10(gags1$P.Value)>=1.3,]$label <- rownames(gags1[abs(gags1$logFC)>=0.3 & -log10(gags1$P.Value)>=1.3,])

res$label <- rownames(res)

ggplot(res, aes(x=logFC,y=-log10(P.Value), label=label)) + 
        geom_point(alpha=0.4, color="lightgray") + theme_bw() +
        geom_point(data=gags1, mapping=aes(x=logFC,y=-log10(P.Value),color=color), size=3) + 
        geom_hline(yintercept = 1.3) +
        geom_vline(xintercept = c(-0.3,0.3)) +
        xlim(c(-1.5, 1.5))+ ylim(c(0,15)) +
        geom_text_repel(data = gags1, label=gags1$label,size = 6, box.padding = 0.5)

##### GENERATE BOX PLOTS used in figures 1, 3, and supplemental figures.


graphs <- list()
#generate one plot per GAG type
for (i in unique(gags$GAG)){
  
  data <- dataset[tolower(rownames(dataset)) %in% tolower(gags[gags$GAG==i,]$Gene),]
  datat <- t(data) %>% as.data.frame() %>% rownames_to_column("sampleID") 
  datat$disease <- factor(datasetmeta$disease_state, levels = c("Normal", "COPD"))
  datat %>%  gather(key="Gene", value = "FPKM",-sampleID, -disease) -> datat1
  
  graphs[[i]] <- ggplot(datat1, aes(x=disease, y=FPKM, fill=disease)) + 
                    facet_wrap(~Gene, scales = "free",ncol=5) + 
                    geom_jitter()+ scale_y_continuous(expand = c(0, 0, .25, 0)) +
                    geom_signif(comparisons = list(c("Normal", "COPD")), step_increase = 2, map_signif_level = T) + 
                    ggtitle(paste0(i, " related genes")) + 
                    geom_boxplot(alpha=0.5) + theme_bw()
  
}
graphs$CSDS
graphs$HS
graphs$CORE

###### Figure 3
plotGeneTPM("CHST11", tpm)

#### Gene Set Enrichment Analysis (Figure 4):
results <- results[order(-results$logFC),]

gene_list <- results$logFC

names(gene_list) <- rownames(results)


library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
gse <- gseGO(gene_list,
             ont = "BP",
             keyType = "SYMBOL",
             OrgDb = "org.Hs.eg.db", 
             eps = 1e-300)

gseaplot(gse, geneSetID = "GO:0007179", title = "GO:0007179: TGF-beta receptor signaling pathway")


###########################################################################################
##### Dataset 2: GSEGSE124180
###########################################################################################

ds2 <- read.csv("GSE124180/GSE124180_gene_count_table.tsv", sep = "\t")
rownames(ds2) <- ds2$ENSEMBL_GENEID
ds2$ENSEMBL_GENEID <- NULL

ds2meta <- read.csv("GSE124180/GSE124180_series_matrix.txt",skip = 27, sep = "\t")
ds2meta <- t(ds2meta)
ds2meta$sampleID <- colnames(ds2)
colnames(ds2meta) <- ds2meta[1,]

ds2meta <- ds2meta[2:nrow(ds2meta),] %>% as.data.frame()
colnames(ds2meta) <- gsub("!", "", colnames(ds2meta))
ds2m <- ds2meta[,c(1,7,11,14,15,16)]
colnames(ds2m) <- c("accession", "type", "age", "packyears", "smoker", "diagnosis")

ds2m$sampleID <- rownames(ds2m)
ds2m <- apply(ds2m, 2, function(x) gsub(".*:","",x)) %>% as.data.frame()

library(annotables)
genes <- rownames(ds2)
ds2$X <- rownames(ds2)
ds2 <- ds2 %>% inner_join(grch38, by = c("X" = "ensgene"))
ds2 <- ds2[2:nrow(ds2),]
ds2 <- na.omit(ds2)
ds2m <- ds2m %>% filter(type=="large airway")

### Generate comparisons with DEseq2
#load libraries 
library(ggrepel)
library(ggbiplot)
library(RColorBrewer)
library(DESeq2)

#load sampleInfo, and data in TPM and read counts.

ds2s <- ds2 %>% select(ds2m$sampleID) 
ds2s <- apply(ds2s, 2, as.numeric)
rownames(ds2s) <- ds2$symbol
dds <- DESeqDataSetFromMatrix(ds2s, colData = ds2m, design = ~diagnosis)
dds <- estimateSizeFactors( dds )
sizeFactors( dds )
dds <-estimateDispersions(dds)
keep <- rowSums(counts(dds))>=1
dds <- dds[keep,]

levels(dds$diagnosis)

dds <- DESeq(dds)

results <- results(dds,alpha=0.99)
summary(results)


# Calculate TPMs

gene_info <- ds2[,64:72]
counts <- ds2[,1:63]
counts <- as.data.frame(apply(counts, 2, as.numeric))

rownames(counts) <- gene_info$entrez

#first get Length in kilobase
geneLength <- (as.numeric(gene_info$end)-as.numeric(gene_info$start))/1000

# divide all counts by the gene Length
TPM <- apply(counts, 2, function(x) {x/geneLength})

# count all RPK sum and divide by 10^6
RPKSUM <- colSums(TPM)/1000000

#Finally, to get TPM divided each RPK by RPKSUM per milion
TPM <- apply(TPM, 1, function(x) {x/RPKSUM})
TPM <- as.data.frame(t(TPM))
TPM$ext_gene <- gene_info$symbol

ds2m[ds2m$diagnosis==" cont",]$diagnosis <- "Normal"
ds2m[ds2m$diagnosis==" case",]$diagnosis <- "COPD"
