# Load libraries
library(DESeq2)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db)
library(biomaRt)
library(BiocManager)
library(AnnotationDbi)
library(clusterProfiler)
library(EnhancedVolcano)
library(dplyr)
library(tibble)
library(ggplot2)
library(enrichplot)

# Set working directory
setwd(
  "/data/users/acastro/breast_cancer/analysis/feature_counts/"
)

##############################   Load Data Frames     ################################

# Load Feature Counts data
counts <-
  read.table("clean_data_counts.txt", header = TRUE, sep = "\t")
head(counts)
dim(counts)

############################## Exploratory Data Analysis #############################

# Data preparation for visualizing gene expression patterns following a multifacor design

# Create a data frame for gene IDs
gen.id <- counts$Geneid
head(gen.id)

# Convert data counts into a matrix
counts <- counts[, 2:13]
rownames(counts) <- gen.id
counts <- as.matrix(counts)
head(counts)
dim(counts) # Visualize 12 samples

# Assign conditions, the experimental groups of cancer cells
coldata <- data.frame("condition" = as.factor(c(
  rep("HER2", 3),
  rep("NonTNBC", 3),
  rep("Normal", 3),
  rep("TNBC", 3)
)),
row.names = colnames(counts))
head(coldata)
dim(coldata)

# Create a DESeqDataSet object
# Normalize to identify the differential expressions genes
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# Data transformation and visualisation
vsd <- vst(dds, blind = TRUE)
vsd
rld <- rlog(dds, blind = TRUE)


# This gives log2(n + 1)
ntd <- normTransform(dds)


################### Data quality assessment by sample clustering and visualization #############

# Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)
head(sampleDistMatrix)

rownames(sampleDistMatrix) <- paste(rld$condition,
                                    rld$type,
                                    sep = "-")

png(
  "heatmap1.png",
  width = 8,
  height = 6,
  units = "in",
  res = 300
) # Export plot
colors <-
  colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors
)
dev.off()


########################## Principal Component Analysis #####################################

# Principal Component Analysis with log2 transformed normalized data
plotPCA <- plotPCA(rld, intgroup = "condition", ntop = 500)

ggsave("PCA_plot.png",
       plot = plotPCA,
       dpi = 300) # Save plot
graphics.off()


########################## Differential expression analysis DESeq ###########################

# Differential expression analysis
dds <- DESeq(dds)

########################## DESeq Contrast between "NonTNBC" and "HER2" ######################

# Results by contrasting conditions "NonTNBC" and "HER2"

# DESeq
res.nonTNBCvsHER <- results(
  dds,
  contrast = c("condition", "NonTNBC", "HER2"),
  alpha = 0.05,
  independentFiltering = TRUE,
  cooksCutoff = TRUE
)

res.nonTNBCvsHER
summary(res.nonTNBCvsHER)

# Filter for genes with nonzero total read count
res.nonTNBCvsHER.f <-
  res.nonTNBCvsHER[res.nonTNBCvsHER$baseMean > 0,]

# Add name column ENSMBL
res.nonTNBCvsHER.f <- as.data.frame(res.nonTNBCvsHER.f) %>%
  rownames_to_column("ENSEMBL")
res.nonTNBCvsHER.f

# Connect to Ensembl and retrieve data with hgnc_symbol
ensembl <- useEnsembl(biomart = "ensembl", host = "asia")

listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
head(datasets)
searchDatasets(mart = ensembl, pattern = "hsapiens")
ensembl <-
  useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
ensembl <-
  useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
listEnsemblArchives()
listEnsembl(version = 110)

ensembl110 <- useEnsembl(biomart = "genes",
                         dataset = "hsapiens_gene_ensembl",
                         version = 110)
ensembl110

gene_data1 <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = ensembl110,
  filters = "ensembl_gene_id",
  values = res.nonTNBCvsHER.f$ENSEMBL
)

gene_data1

# Convert ensembl_gene_id to character in gene_data
gene_data1 <-
  gene_data1 %>% mutate(ensembl_gene_id = as.character(ensembl_gene_id))
gene_data1

# Join results with gene_data
sig.nonTNBCvsHER.j <- left_join(res.nonTNBCvsHER.f,
                                gene_data1,
                                by = c("ENSEMBL" = "ensembl_gene_id"))

sig.nonTNBCvsHER.s <-
  sig.nonTNBCvsHER.j[order(sig.nonTNBCvsHER.j$padj),]
sig.nonTNBCvsHER.s

sig.nonTNBCvsHER.s
write.csv(sig.nonTNBCvsHER.j, "sig.NonTNBCvsHER2.j.csv")
dim(sig.nonTNBCvsHER.s)

# Visualization of DEGs
volcano_plot_NonTNBCvsHER <- EnhancedVolcano(
  sig.nonTNBCvsHER.s,
  lab = sig.nonTNBCvsHER.s$hgnc_symbol,
  selectLab = sig.nonTNBCvsHER.s[c(1:15), 8],
  x = "log2FoldChange",
  y = "pvalue",
  title = "NonTNBC vs HER2",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.5,
  labSize = 6.0,
  xlim = c(-20, 20)
)
ggsave("volcano_plot_NonTNBCvsHER.png",
       plot = volcano_plot_NonTNBCvsHER,
       dpi = 300) # Save plot
graphics.off()


# Genes are differential expressed (DE) in the pairwise comparison padj < 0.5
sig.nonTNBCvsHER <- filter(sig.nonTNBCvsHER.s, padj < 0.05)
write.csv(sig.nonTNBCvsHER, "sig.NonTNBCvsHER2.j.csv")

sig.nonTNBCvsHER
dim(sig.nonTNBCvsHER)


########################## DESeq Contrast between "TNBC" and "HER2" ########################

# Results by contrasting conditions "TNBC" and "HER2"

res.TNBCvsHER <- results(
  dds,
  contrast = c("condition", "TNBC", "HER2"),
  alpha = 0.05,
  independentFiltering = TRUE,
  cooksCutoff = TRUE
)

summary(res.TNBCvsHER)

# Filter for genes with nonzero total read count
res.TNBCvsHER.f <- res.TNBCvsHER[res.TNBCvsHER$baseMean > 0,]
dim(res.TNBCvsHER.f)

# Add name column ENSMBL
res.TNBCvsHER.f <- as.data.frame(res.TNBCvsHER.f) %>%
  rownames_to_column("ENSEMBL")
res.TNBCvsHER.f
dim(res.TNBCvsHER.f)

# Connect to Ensembl and retrieve data with hgnc_symbol
gene_data2 <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = ensembl110,
  filters = "ensembl_gene_id",
  values = res.TNBCvsHER.f$ENSEMBL
)


# Convert ensembl_gene_id to character in gene_data
gene_data2 <-
  gene_data2 %>% mutate(ensembl_gene_id = as.character(ensembl_gene_id))
gene_data2
dim(gene_data2)


# Join results with gene_data
sig.TNBCvsHER.j <- left_join(res.TNBCvsHER.f,
                             gene_data2,
                             by = c("ENSEMBL" = "ensembl_gene_id"))
dim(sig.TNBCvsHER.j)

# order padj value
sig.TNBCvsHER.s <- sig.TNBCvsHER.j[order(sig.TNBCvsHER.j$padj),]
sig.TNBCvsHER.s
dim(sig.TNBCvsHER.s)

# Visualization of DEGs
volcano_plot_TNBCvsHER <- EnhancedVolcano(
  sig.TNBCvsHER.s,
  lab = sig.TNBCvsHER.s$hgnc_symbol,
  selectLab = sig.TNBCvsHER.s[c(1:45), 8],
  x = "log2FoldChange",
  y = "pvalue",
  title = "TNBC vs HER2",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.5,
  labSize = 6.0,
  xlim = c(-20, 20)
)
ggsave("volcano_plot_TNBCvsHER2.png",
       plot = volcano_plot_TNBCvsHER,
       dpi = 300) # Save plot
graphics.off()

# Genes are differential expressed (DE) in the pairwise comparison padj < 0.5
sig.TNBCvsHER <- filter(sig.TNBCvsHER.s, padj < 0.05)
sig.TNBCvsHER
write.csv(sig.TNBCvsHER, "sig.TNBCvsHER2.j.csv")
dim(sig.TNBCvsHER)

########################## DESeq Contrast between "NonTNBC" and "TNBC"  ###################

res.NonTNBCvsTNBC <- results(
  dds,
  contrast = c("condition", "NonTNBC", "TNBC"),
  alpha = 0.05,
  independentFiltering = TRUE,
  cooksCutoff = TRUE
)

summary(res.NonTNBCvsTNBC)

# Filter for genes with nonzero total read count
res.NonTNBCvsTNBC.f <-
  res.NonTNBCvsTNBC[res.NonTNBCvsTNBC$baseMean > 0,]
dim(res.NonTNBCvsTNBC.f)

# Add name column ENSMBL
res.NonTNBCvsTNBC.f <- as.data.frame(res.NonTNBCvsTNBC.f) %>%
  rownames_to_column("ENSEMBL")
res.NonTNBCvsTNBC.f
dim(res.NonTNBCvsTNBC.f)

# Connect to Ensembl and retrieve data with hgnc_symbol
gene_data2 <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = ensembl110,
  filters = "ensembl_gene_id",
  values = res.NonTNBCvsTNBC.f$ENSEMBL
)


# Convert ensembl_gene_id to character in gene_data
gene_data2 <-
  gene_data2 %>% mutate(ensembl_gene_id = as.character(ensembl_gene_id))
gene_data2
dim(gene_data2)


# Join results with gene_data
sig.NonTNBCvsTNBC.j <- left_join(res.NonTNBCvsTNBC.f,
                                 gene_data2,
                                 by = c("ENSEMBL" = "ensembl_gene_id"))
dim(sig.NonTNBCvsTNBC.j)

# order padj value
sig.NonTNBCvsTNBC.s <-
  sig.NonTNBCvsTNBC.j[order(sig.NonTNBCvsTNBC.j$padj),]
sig.NonTNBCvsTNBC.s
dim(sig.NonTNBCvsTNBC.s)

# Visualization of DEGs
volcano_plot_NonTNBCvsTNBC <- EnhancedVolcano(
  sig.NonTNBCvsTNBC.s,
  lab = sig.NonTNBCvsTNBC.s$hgnc_symbol,
  selectLab = sig.NonTNBCvsTNBC.s[c(1:15), 8],
  x = "log2FoldChange",
  y = "pvalue",
  title = "NonTNBC vs TNBC",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.5,
  labSize = 6.0,
  xlim = c(-20, 20)
)
ggsave("volcano_plot_NonTNBCvsTNBC.png",
       plot = volcano_plot_NonTNBCvsTNBC,
       dpi = 300) # Save plot
graphics.off()

# Genes are differential expressed (DE) in the pairwise comparison padj < 0.5
sig.NonTNBCvsTNBC <- filter(sig.NonTNBCvsTNBC.s, padj < 0.05)
sig.NonTNBCvsTNBC
write.csv(sig.NonTNBCvsTNBC, "sig.NonTNBCvsTNBC.j.csv")
dim(sig.NonTNBCvsTNBC)


########################## DESeq Contrast between "TNBC" and "Normal"  ###################

res.TNBCvsNormal <- results(
  dds,
  contrast = c("condition", "TNBC", "Normal"),
  alpha = 0.05,
  independentFiltering = TRUE,
  cooksCutoff = TRUE
)

summary(res.TNBCvsNormal)

# Filter for genes with nonzero total read count
res.TNBCvsNormal.f <-
  res.TNBCvsNormal[res.TNBCvsNormal$baseMean > 0,]
dim(res.TNBCvsNormal.f)

# Add name column ENSMBL
res.TNBCvsNormal.f <- as.data.frame(res.TNBCvsNormal.f) %>%
  rownames_to_column("ENSEMBL")
res.TNBCvsNormal.f
dim(res.TNBCvsNormal.f)

# Connect to Ensembl and retrieve data with hgnc_symbol
gene_data2 <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = ensembl110,
  filters = "ensembl_gene_id",
  values = res.TNBCvsNormal.f$ENSEMBL
)


# Convert ensembl_gene_id to character in gene_data
gene_data2 <-
  gene_data2 %>% mutate(ensembl_gene_id = as.character(ensembl_gene_id))
gene_data2
dim(gene_data2)


# Join results with gene_data
sig.TNBCvsNormal.j <- left_join(res.TNBCvsNormal.f,
                                gene_data2,
                                by = c("ENSEMBL" = "ensembl_gene_id"))
dim(sig.TNBCvsNormal.j)

# order padj value
sig.TNBCvsNormal.s <-
  sig.TNBCvsNormal.j[order(sig.TNBCvsNormal.j$padj),]
sig.TNBCvsNormal.s
dim(sig.TNBCvsNormal.s)

# Visualization of DEGs
volcano_plot_TNBCvsNormal <- EnhancedVolcano(
  sig.TNBCvsNormal.s,
  lab = sig.TNBCvsNormal.s$hgnc_symbol,
  selectLab = sig.TNBCvsNormal.s[c(1:15), 8],
  x = "log2FoldChange",
  y = "pvalue",
  title = "TNBC vs Normal",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.5,
  labSize = 6.0,
  xlim = c(-20, 20)
)
ggsave("volcano_plot_TNBCvsNormal.png",
       plot = volcano_plot_TNBCvsNormal,
       dpi = 300) # Save plot
graphics.off()

# Genes are differential expressed (DE) in the pairwise comparison padj < 0.5
sig.TNBCvsNormal <- filter(sig.TNBCvsNormal.s, padj < 0.05)
sig.TNBCvsNormal
write.csv(sig.TNBCvsNormal, "sig.TNBCvsNormal.j.csv")
dim(sig.TNBCvsNormal)


########################## DESeq analysis genes selected ########################

# Selected genes were evaluated to analyze number of normalized counts between contrasts. 
# List of genes to analyze
genes_to_analyze <- c("RIMS4", "CDKN2A", "KCNK15", "HRCT1", "FSIP1", "FOXA1")

# Use normalized dataset
normalized_counts <- counts(dds, normalized = TRUE)

# Make matrix into a dataframe
normalized_counts_df <- as.data.frame(normalized_counts) %>%
  rownames_to_column("ENSEMBL")

# Connect to Ensembl and retrieve data with hgnc_symbol
gene_data.norm <-
  getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    mart = ensembl110,
    filters = "ensembl_gene_id",
    values = normalized_counts_df$ENSEMBL
  )

# Convert ensembl_gene_id to character in gene_data
gene_data.norm <-
  gene_data.norm %>% mutate(ensembl_gene_id = as.character(ensembl_gene_id))

# Join results with gene_data
merged_data.norm <-
  left_join(normalized_counts_df,
            gene_data.norm,
            by = c("ENSEMBL" = "ensembl_gene_id"))

# Function to create bar plots
create_gene_bar_plot <- function(gene_name) {
  # Extract gene expression data for the specific gene
  gene_index <- which(merged_data.norm$hgnc_symbol == gene_name)
  gene_data <- merged_data.norm[gene_index, grep("HER|NonTNBC|TNBC|Normal",
                                                 colnames(merged_data.norm),
                                                 value = TRUE)]
  
  plot_data <- data.frame(Sample = names(gene_data),
                          Expression = as.numeric(gene_data),
                          Gene = gene_name)
  
  # Create the bar plot
  gg_plot <- ggplot(plot_data, aes(x = Sample, y = Expression, fill = Sample)) +
    geom_bar(stat = "identity") +
    labs(title = gene_name, y = "Normalized counts") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  return(gg_plot)
}

# Create a list to store individual plots
gene_plots <- list()

# Loop through selected genes
for (gene_name in genes_to_analyze) {
  gene_plots[[gene_name]] <- create_gene_bar_plot(gene_name)
}

# Combine plots into a single figure
combined_plot <- cowplot::plot_grid(plotlist = gene_plots, ncol = 2, byrow = TRUE)

# Save the combined plot
ggsave(
  "combined_gene_plots.png",
  combined_plot,
  width = 9,
  height = 12,
  dpi = 300
)


########################## Overrepresentation analysis "NonTNBC" and "HER2"##################

## Create background dataset for hypergeometric testing using all genes tested for significance in the results

# Data of differential expression total
sig.nonTNBCvsHER.s

# Extract ENSEBL column
allOEgenes_1 <- as.character(sig.nonTNBCvsHER.s$ENSEMBL)
allOEgenes_1
write.csv(allOEgenes_1, "allOEgenes_NonTNBCvsHER.csv")

## Significant results of DESeq
sigOEgenes_1 <- as.character(sig.nonTNBCvsHER$ENSEMBL)
sigOEgenes_1
write.csv(sigOEgenes_1, "sigOEgenes_NonTNBCvsHER.csv")

## Run GO enrichment analysis
ego1 <- enrichGO(
  gene = sigOEgenes_1,
  universe = allOEgenes_1,
  keyType = "ENSEMBL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE
)

## Output results from GO analysis to a table
cluster_summary1 <- data.frame(ego1)
write.csv(cluster_summary1, "cluster_summaryNonTNBCvsHER.csv")

# Visualizing clusterProfiler results
# Dotplot
dotplot(ego1, showCategory = 8) + ggtitle("NonTNBC vs HER2")
ggsave(
  "dotploNonTNBCvsHER2.png",
  plot = last_plot(),
  device = "png",
  width = 7,
  height = 10
)
graphics.off()

# bar plot
barplot(ego1, showCategory = 8)
ggsave(
  "barplotenonTNBCvsHER2.png",
  plot = last_plot(),
  device = "png",
  width = 7,
  height = 10
)
graphics.off()

# bar plot
mutate(ego1, qscore = -log(p.adjust, base = 10)) %>%
  barplot(x = "qscore")
ggsave(
  "qscoreNonTNBCvsHER2.png",
  plot = last_plot(),
  device = "png",
  width = 7,
  height = 10
)
graphics.off()

# Gene-Concept Network
edox1 <- setReadable(ego1, "org.Hs.eg.db", "ENSEMBL")
edox1

png(
  "cnetplononTNBCvsHER2.png",
  width = 4000,
  height = 3500,
  res = 300
)
cnetplot(edox1, circular = TRUE, colorEdge = TRUE)
dev.off()

########################## Overrepresentation analysis "TNBC" and "HER2"#####################

## Create background dataset for hypergeometric testing using all genes tested for significance in the results

# Data of differential expression total sig.TNBCvsHER.s
sig.TNBCvsHER.s

# Extract ENSEBL column
allOEgenes_2 <- as.character(sig.TNBCvsHER.s$ENSEMBL)
allOEgenes_2
write.csv(allOEgenes_2, "allOEgenes_TNBCvsHER.csv")


## Significant results of DESeq
sigOEgenes_2 <- as.character(sig.TNBCvsHER$ENSEMBL)
sigOEgenes_2
write.csv(sigOEgenes_2, "sigOEgenes_TNBCvsHER.csv")


## Run GO enrichment analysis
ego2 <- enrichGO(
  gene = sigOEgenes_2,
  universe = allOEgenes_2,
  keyType = "ENSEMBL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE
)

## Output results from GO analysis to a table
cluster_summary2 <- data.frame(ego2)
write.csv(cluster_summary2, "cluster_summaryTNBCvsHER.csv")

# Visualizing clusterProfiler results
# Dotplot
dotplot(ego2, showCategory = 10) + ggtitle("TNBC vs HER2")
ggsave(
  "dotplotTNBCvsHER2.png",
  plot = last_plot(),
  device = "png",
  width = 7,
  height = 10
)
graphics.off()

# bar plot
barplot(ego2, showCategory = 10)
ggsave(
  "barplotTNBCvsHER2.png",
  plot = last_plot(),
  device = "png",
  width = 7,
  height = 10
)
graphics.off()

# bar plot
mutate(ego2, qscore = -log(p.adjust, base = 10)) %>%
  barplot(x = "qscore")
ggsave(
  "qscoreTNBCvsHER2.png",
  plot = last_plot(),
  device = "png",
  width = 7,
  height = 10
)
graphics.off()

# Gene-Concept Network
edox2 <- setReadable(ego2, "org.Hs.eg.db", "ENSEMBL")
edox2

png(
  "cnetplotTNBCvsHER2.png",
  width = 4000,
  height = 3500,
  res = 300
)
cnetplot(edox2, circular = TRUE, colorEdge = TRUE)
dev.off()

########################## Overrepresentation analysis "NonTNBC" and "TNBC"#####################

# Data of differential expression total sig.NonTNBCvsHER.s
sig.NonTNBCvsTNBC.s

# Extract ENSEBL column
allOEgenes_3 <- as.character(sig.NonTNBCvsTNBC.s$ENSEMBL)
allOEgenes_3
write.csv(allOEgenes_3, "allOEgenes_NonTNBCvsTNBC.csv")


## Significant results of DESeq
sigOEgenes_3 <- as.character(sig.NonTNBCvsTNBC$ENSEMBL)
sigOEgenes_3
write.csv(sigOEgenes_3, "sigOEgenes_NonTNBCvsTNBC.csv")


## Run GO enrichment analysis
ego3 <- enrichGO(
  gene = sigOEgenes_3,
  universe = allOEgenes_3,
  keyType = "ENSEMBL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE
)

## Output results from GO analysis to a table
cluster_summary3 <- data.frame(ego3)
write.csv(cluster_summary3, "cluster_summaryNonTNBCvsTNBC.csv")

# Visualizing clusterProfiler results
# Dotplot
dotplot(ego3, showCategory = 10) + ggtitle("NonTNBC vs TNBC")
ggsave(
  "dotplotNonTNBCvsTNBC.png",
  plot = plot3,
  device = "png",
  width = 7,
  height = 15
)
graphics.off()

# bar plot
barplot(ego3, showCategory = 10)
ggsave(
  "barplotNonTNBCvsTNBC.png",
  plot = last_plot(),
  device = "png",
  width = 7,
  height = 10
)
graphics.off()

# bar plot
mutate(ego3, qscore = -log(p.adjust, base = 10)) %>%
  barplot(x = "qscore")
ggsave(
  "qscoreNonTNBCvsTNBC.png",
  plot = last_plot(),
  device = "png",
  width = 7,
  height = 10
)
graphics.off()

# Gene-Concept Network
edox3 <- setReadable(ego3, "org.Hs.eg.db", "ENSEMBL")
edox3

png(
  "cnetplotNonTNBCvsTNBC.png",
  width = 4000,
  height = 3500,
  res = 300
)
cnetplot(edox3, circular = TRUE, colorEdge = TRUE)
dev.off()

# Generate a single figure combining plots 
plot1 <- dotplot(ego1, showCategory = 8) + ggtitle("NonTNBC vs HER2")
plot2 <- dotplot(ego2, showCategory = 10) + ggtitle("TNBC vs HER2")
plot3 <- dotplot(ego3, showCategory = 10) + ggtitle("NonTNBC vs TNBC")

cnet_plot1 <- cnetplot(edox1, circular = TRUE, colorEdge = TRUE)
cnet_plot2 <- cnetplot(edox2, circular = TRUE, colorEdge = TRUE)
cnet_plot3 <- cnetplot(edox3, circular = TRUE, colorEdge = TRUE)

# Combine dotplot and cnet plots
combined_plots <- (plot1 | cnet_plot1) / (plot2 | cnet_plot2) / (plot3 | cnet_plot3)

# Save the combined plot
ggsave("combined_plots.png", plot = combined_plots, width = 30, height = 30, dpi = 300, units = "in")
graphics.off()


######### Bar plot from the output results of summarized feature counts ############

# Create a data frame with the data
data <- data.frame(
  Status = c("Assigned", "Unassigned_Unmapped", "Unassigned_MultiMapping", "Unassigned_NoFeatures", "Unassigned_Ambiguity"),
  HER21 = c(10.9, 2.3, 143.1, 15.0, 1.5),
  HER22 = c(15.6, 1.6, 146.5, 17.6, 1.8),
  HER23 = c(11.0, 0.9, 121.9, 11.8, 1.6),
  NonTNBC1 = c(21.0, 2.7, 122.0, 9.7, 3.2),
  NonTNBC2 = c(11.8, 2.3, 114.4, 10.3, 1.7),
  NonTNBC3 = c(12.7, 2.5, 122.5, 11.4, 1.8),
  Normal1 = c(12.4, 0.2, 4.4, 0.8, 1.1),
  Normal2 = c(26.5, 0.4, 6.8, 2.1, 2.0),
  Normal3 = c(31.6, 0.3, 5.3, 1.4, 2.2),
  TNBC1 = c(11.2, 1.6, 74.9, 13.9, 1.3),
  TNBC2 = c(9.3, 2.2, 116.2, 7.7, 1.2),
  TNBC3 = c(10.0, 2.7, 120.4, 8.0, 1.2)
)

data_long <- tidyr::gather(data, key = "Sample", value = "Value", -Status)

# Stacked bar plot 
ggplot(data_long, aes(x = Sample, y = Value, fill = Status)) +
  geom_bar(stat = "identity") +
  labs(title = "B",
       y = "Count (in million)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Assigned" = "#0b7a75", 
                               "Unassigned_Unmapped" = "red", 
                               "Unassigned_MultiMapping" = "#ab4a43", 
                               "Unassigned_NoFeatures" = "#aba194", 
                               "Unassigned_Ambiguity" = "#d39f94"))


ggsave("stacked_bar_plot.png", 
       plot = last_plot() +
         theme_minimal() +  
         theme(panel.background = element_rect(fill = "white"),
               axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),  
               axis.text.y = element_text(size = 13, color = "black"),  
               axis.title = element_text(size = 13, color = "black"),  
               legend.position = "right", 
               legend.text = element_text(size = 13, color = "black"),  
               plot.title = element_text(size = 13, color = "black")  
         ), 
       width = 12,  
       height = 6, 
       dpi = 300,
       bg = "white"  
)


