library(tidyverse)
library(DESeq2)

#### Importing files ####
# Importing count data
countData <- read.csv("counts.csv", header = TRUE, row.names = 1) %>% 
  as.matrix()

# Make changes to transcript names
rownames(countData) <- gsub("Glyma.", "GLYMA_", rownames(countData))

# Filter data where you only have 0 or 1 read count across all samples.
countData = countData[rowSums(countData)>1, ]
glimpse(countData)

# Import metadata
colData = read.csv("meta.csv", row.names=1)
colData
colData <- colData[,-c(1,4)]
# colnames(colData)[2] <- "treatment"
colData$cultivar.treatment <- gsub("H","T",colData$cultivar.treatment)
rownames(colData) <- paste(as.character(colData$cultivar.treatment),
                           c(1:length(colData$cultivar.treatment)),
                           sep = "_")
colData <- colData[,-4]
colData$treatment <- gsub(" ", "_", colData$treatment)

head(countData)
colnames(countData) <- rownames(colData)

# Convert the treatment variables to factors
for (i in 1:length(colnames(colData))){
  colData[,i] <- factor(colData[,i])
}
glimpse(colData)

# Checking if column names of the gene counts match the row names of the metadata
check <- colnames(countData) == row.names(colData)
# Create the DESeqDataSet object if all checks pass
# ?DESeqDataSetFromMatrix()

#### DEG Analysis ####
# Set up the DESeqDataSet Object and run the DESeq pipeline
if (all(check)){
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ cultivar + treatment + cultivar:treatment)
} else{ 
  stop("The column names of the gene counts does not match the row names of the metadata")
}

dds$cultivar
dds$treatment
dds$treatment <- relevel(dds$treatment, "normal_watered")

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~1)
levels(dds$cultivar)
levels(dds$treatment)
dds$treatment <- relevel(dds$treatment, ref = "NW")
dds$cultivar <- relevel(dds$cultivar, ref = "tolerant")
dds$cultivar <- factor(dds$cultivar, levels = c("sensitive", "tolerant"))

mm <- model.matrix(~ treatment + cultivar + treatment:cultivar, colData(dds))
dds <- DESeq(dds, full=mm)

resultsNames(dds)

# Normalizes the count data (median of ratios method of normalization)
# Estimates library size and dispersion
dds = DESeq(dds)
dds
resultsNames(dds)
results(dds, contrast = c("cultivar_cultivarsensitive.treatment_treatmentLD"))

colData(dds)

# Plot the dispersion of the gene estimates
plotDispEsts(dds, main="DESeq: Per-gene dispersion estimates")

# Principal components biplot on variance stabilized data
print(plotPCA(varianceStabilizingTransformation(dds), 
              intgroup=c("treatment", "cultivar")))

print(plotPCA(varianceStabilizingTransformation(dds), 
              intgroup=c("treatment_name")))

print(plotPCA(varianceStabilizingTransformation(dds), 
              intgroup=c("cultivar")))

resultsNames(dds)

results(dds, name = "Intercept")
# Function to extract the contrast you want
extract.contrast <- function(data, column, control = "", treatment){
  dat <- results(data, contrast = 
                   c(column, 
                     control, 
                     treatment)); dat <- dat[order(
                       dat$pvalue),]; summary(dat)
  return(dat)
}
H_SD <- results(dds, contrast = list("cultivartolerant.treatmentNW", 
                                     "cultivartolerant.treatmentSD"))

H_SD <- results(dds, contrast = list("cultivartolerant.treatmentSD",
                                     "cultivartolerant.treatmentNW"))

results(dds, contrast = c("cultivar", "sensitive", "tolerant"))
results(dds, contrast = c("cultivar", "tolerant", "sensitive"))

results(dds, contrast = c("cultivar:treatment", "sensitive", "LD"))

# Extract the contrast you want
S_LD <- extract.contrast(dds, "cultivar.treatment", "S_NW", "S_LD")
S_SD <- extract.contrast(dds, "cultivar.treatment", "S_NW", "S_SD")
H_LD <- extract.contrast(dds, "cultivar.treatment", "S_NW", "H_LD")
H_SD <- extract.contrast(dds, "cultivar.treatment", "S_NW", "H_SD")

#### Adding gene annotations ####
# Getting entrez ids
library("AnnotationDbi")
library(biomaRt)
org.Gmx.eg.db <- loadDb("Gm.db")
columns(org.Gmx.eg.db)

marty <- useMart(host = "https://jul2018-plants.ensembl.org", 
                    biomart = "plants_mart", 
                    dataset = "gmax_eg_gene")

listAttributes(mart)
univGmx.geneID <- getBM(attributes = c("ensembl_gene_id", "entrezgene", "kegg"),mart = marty)

univGmx.geneID <- univGmx.geneID[!is.na(univGmx.geneID$entrezgene),]
univGmx.geneID <- univGmx.geneID[!duplicated(univGmx.geneID$entrezgene),]

H_LD$ENTREZID <- univGmx.geneID$entrezgene[match(
  rownames(H_LD), univGmx.geneID$ensembl_gene_id)]

# Adding entrez ids to the results of the DESeq analysis
add.entrez <- function(target.dat, universal.dat){
  target.dat$ENTREZID <- universal.dat$entrezgene[match(
    rownames(target.dat), universal.dat$ensembl_gene_id)]
  return(target.dat)
}

H_LD <- add.entrez(H_LD, univGmx.geneID)
H_SD <- add.entrez(H_SD, univGmx.geneID)
S_LD <- add.entrez(S_LD, univGmx.geneID)
S_SD <- add.entrez(S_SD, univGmx.geneID)

# Adding new columns to the results of the DESeq analysis
# Function to add annotation columns
add.annotations <- function(orgDB, data){
  data$symbol = mapIds(orgDB,
                    keys = as.character(data$ENTREZID), 
                    column = "SYMBOL",
                    keytype = "ENTREZID",
                    multiVals = "first")
  data$name =   mapIds(orgDB,
                    keys = as.character(data$ENTREZID), 
                    column = "GENENAME",
                    keytype = "ENTREZID",
                    multiVals = "first")
  return(data)
}

# select(org.Gmx.eg.db,column = "ENTREZID")
# 
# keys(org.Gmx.eg.db, "ENTREZID")

# subs <- data.frame(ENTREZID = keys(org.Gmx.eg.db, "ENTREZID"),
#            SYMBOL = keys(org.Gmx.eg.db, "SYMBOL"),
#            GENENAME = keys(org.Gmx.eg.db, "GENENAME"))

S_LD <- add.annotations(org.Gmx.eg.db, S_LD)
S_SD <- add.annotations(org.Gmx.eg.db, S_SD)
H_LD <- add.annotations(org.Gmx.eg.db, H_LD)
H_SD <- add.annotations(org.Gmx.eg.db, H_SD)

glimpse(dat)

dat1 <- dat %>% 
  mutate(symbol = as.character(symbol),
         )

dataSheets <- c("S_LD", "S_SD", "H_LD", "H_SD")

for (i in 1:length(dataSheets)){
  dat <- data.frame(get(dataSheets[i]))
  file <- sprintf("%s.csv", dataSheets[i])
  write.csv(dat, file)
}

class(dat)

library(writexl)
write_xlsx(dat, path = "my_data.xlsx")

#### Pathway analysis ####
# KEGG pathway
library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs) # loading kegg datasets for Homo sapiens
data(sigmet.idx.hs) # loading signaling and metabolic pathways
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs] # subsetting just sig & met
head(kegg.sets.hs, 3)

foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)

# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)

# Look at both up (greater), down (less), and statistics.
lapply(keggres, head)


# Get the top 5 upregulated pathways
keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  as_tibble() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()

keggrespathways

# Get the IDs.
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

# Plotting pathways using pathview
# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))


### Gene Ontology ###
data(go.sets.hs) # Loading all GO terms for Homo sapiens
data(go.subs.hs) # Loading all indexes for BP, MF, and CC for Homo sapiens
gobpsets = go.sets.hs[go.subs.hs$BP] # subsetting just BP

# Get the results
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)

# Look at the results for both up (greater), down (less), and statistics.
lapply(gobpres, head)


#### Test trials ####

# data$entrez = mapIds(orgDB,
#                      keys=row.names(data), 
#                      column="ENTREZID",
#                      keytype="ENSEMBL",
#                      multiVals="first")

# H_LD$symbol <- mapIds(org.Gmx.eg.db, keys = as.character(H_LD$ENTREZID),
#                       column = "SYMBOL",
#                       keytype = "ENTREZID",
#                       multiVals = "first")
# 
# keytypes(org.Gmx.eg.db)
# head(keys(org.Gmx.eg.db, "PMID"))