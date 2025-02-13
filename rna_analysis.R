# rna analysis

library(org.Hs.eg.db)
library(clusterProfiler)
library(GSEABase)
library(GSVA)
library(edgeR)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(dplyr)
library(ComplexHeatmap)
library(limma)
library(pheatmap)
library(RColorBrewer)
library(enrichplot)


# RNA expression correlation

##  TMM normalize
anno = read.delim("anno.txt", row.names = 1)

### load data
counts <- read.delim("RNA_readcount.txt", sep="\t")
colnames(counts)[3:14] <- colnames(counts)[3:14] %>% gsub("\\.","-",.)

counts_dedup <- counts[,c(2:14)]
counts_dedup <- counts_dedup[!duplicated(counts_dedup$Gene_id), ]

### normalize
group=anno$Type
count_input <- counts_dedup[,c(anno$RNA_ID)]
row.names(count_input) <- counts_dedup$Gene_id

y <- DGEList(counts=count_input, group=group)

keep <- filterByExpr(y, group=group)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
bcv <- 0.4  # well-controlled experiments are 0.4 for human data

et <- exactTest(y, dispersion=bcv^2)
top <- topTags(et)

tmm_log <- cpm(y, log=TRUE)

exp.df <- t(tmm_log)
exp.df <- merge(anno, exp.df, by="row.names")

## draw
pdf("RNA_cor_per_patient.pdf", 16, 4)
exp.df %>% filter(Sample %in% c("GBO-001", "GBO-002", "GBO-003", "GBO-005", "GBO-008")) %>% column_to_rownames("RNA_ID") %>%
  select(Type, Patient, colnames(exp.df)[23:21125]) %>% pivot_longer(colnames(exp.df)[23:21125], names_to = "Gene", values_to = "Expression") %>%
  pivot_wider(names_from = Type , values_from = Expression ) %>% 
  ggplot(aes(x=Tissue, y=Organoid)) +  
  stat_density_2d(aes(fill = stat(level)), alpha=0.9, geom = "polygon", show.legend = T , color="black", size = 0.2, bins = 9) + 
  scale_fill_distiller(palette = "Spectral") +
  theme_base() + stat_cor(method = "spearman") +
  facet_wrap(~ Patient, scales = "free", ncol = 5) + labs(x = "mRNA expression (Tissue)", y = "mRNA expression (Organoid)")
dev.off()


# KEGG pathway GSVA

## load data
tpm_final <- read.delim("RNA_TPM.txt")
colnames(tpm_final)[3:15] <- colnames(tpm_final)[3:15] %>% gsub("\\.","-",.)
tpm_final.df <- tpm_final[,c(2:15)]
tpm_final.df <- tpm_final.df[!duplicated(tpm_final.df$Gene_id), ]

tpm_final_mat = as.matrix(tpm_final.df %>% remove_rownames() %>% column_to_rownames("Gene_id"))

kegg_geneSets = getGmt("c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt")

## run GSVA
tpm_final_mat=as.matrix(tpm_final.df %>% remove_rownames() %>% column_to_rownames("Gene_id"))

kegg_gsva=gsva(tpm_final_mat, kegg_geneSets, method="gsva",  mx.diff=TRUE, verbose=FALSE, parallel.sz=1)

kegg_gsva=data.frame(t(kegg_gsva), check.names = F, stringsAsFactors = F)

gbo_rna_names = c("GBO-001-OR3-R", "GBO-002-OR5-R", "GBO-003-OR4-R", "GBO-005-OR4-R", "GBO-006-OR7-R", "GBO-008-OR-p-3-R", "GBO-019-OR2-R", "GBO-015-OR-R")

## draw

pdf("KEGG_MMR_barplot_8GBO.pdf", 5, 5)

kegg_gsva[gbo_rna_names,c("KEGG_MEDICUS_REFERENCE_MISMATCH_REPAIR", "KEGG_MEDICUS_REFERENCE_HOMOLOGOUS_RECOMBINATION")] %>%
  rownames_to_column("Sample") %>%
  pivot_longer(c("KEGG_MEDICUS_REFERENCE_MISMATCH_REPAIR", "KEGG_MEDICUS_REFERENCE_HOMOLOGOUS_RECOMBINATION"),
               names_to = "KEGG", values_to = "GSVA") %>%
  mutate(KEGG = gsub("KEGG_MEDICUS_REFERENCE_MISMATCH_REPAIR","Mismatch Repair", KEGG)) %>%
  filter(KEGG == "Mismatch Repair") %>% 
  mutate(Sample = substr(Sample, 1, 7)) %>%
  ggplot(aes(x = Sample, y = GSVA)) + 
  geom_bar(position=position_dodge(), stat="identity", color="black") + 
  theme_base(base_size = 15) +
  theme(axis.text.x=element_text(color="black",angle = 45, hjust=1))+
  labs(x = NULL, y = "GSVA score", subtitle = "KEGG Mismatch Repair") + 
  geom_hline(yintercept=0, size=0.1) + 
  theme(panel.grid.major.y = element_line(color = 'gray80', linewidth = 0.5))

dev.off()


# enrichment analysis

## filter gene
tpm_final.df %>% select("Gene_id", "GBO-005-OR4-R", "GBO-019-OR2-R") %>% remove_rownames() %>% 
  column_to_rownames("Gene_id") %>%
  apply(2, function(x) {x > 5}) %>% as.data.frame %>% rowSums() %>% 
  as.data.frame %>% rename(exp_thres = ".") %>% 
  filter(exp_thres > 0) %>% row.names() -> gbo005_over5_gene_id

## select gene
tpm_final.df %>% select("Gene_id", "GBO-005-OR4-R", "GBO-019-OR2-R") %>% filter(Gene_id %in% gbo005_over5_gene_id) %>%
  mutate(GBO005 = log2(`GBO-005-OR4-R` + 1)) %>% mutate(GBO005R = log2(`GBO-019-OR2-R` + 1)) %>%
  mutate(FC = GBO005R - GBO005) %>% filter(FC >= 2) %>% arrange(desc(FC)) %>% select(Gene_id) %>% deframe -> GBO5R_fc_over2

## enrichment analysis
GBO5R_fc_over2_ego <-enrichGO(GBO5R_fc_over2, OrgDb='org.Hs.eg.db', keyType = "SYMBOL", ont = "ALL",
                         pvalueCutoff = 0.05, pAdjustMethod = "BH",
                         qvalueCutoff = 0.25, minGSSize = 10, maxGSSize = 500)

## draw
pdf("GBO-019_FC_ego.dotplot.pdf", 7.5, 7)
dotplot(GBO5R_fc_over2_ego %>% filter(ONTOLOGY == "BP"), showCategory = 15, label_format = 50, 
        title="GBO-019, Biologial process")
dev.off()


# MGMT mRNA expression

tpm.anno.df <- merge(anno, tpm_final.df %>% remove_rownames() %>% column_to_rownames("Gene_id") %>% t() , by="row.names")

pdf("MGMT_mRNA_barplot.pdf", 4.5, 4)
tpm.anno.df %>% filter(Type == "Organoid") %>% remove_rownames() %>% column_to_rownames("Row.names") %>%
  select("Sample", "MGMT") %>% 
  ggplot(aes(x = Sample, y = MGMT)) + 
  geom_bar(position=position_dodge(), stat="identity", color="black") + 
  theme_base(base_size = 15) +
  theme(axis.text.x=element_text(color="black",angle = 45, hjust=1))+
  geom_hline(yintercept=0, size=0.1) + 
  labs(x = NULL, y = "TPM", subtitle = "MGMT RNA epxression") +
  theme(panel.grid.major.y = element_line(color = 'gray80', linewidth = 0.5))
dev.off()


# VEGFA, VEGFR expression

pdf("VEGFA_VEGFR_tpm.pdf", 7, 4)
tpm_final.df %>% filter(Gene_id %in% c("VEGFA", "FLT1", "KDR", "FLT4")) %>% 
  pivot_longer(colnames(tpm_final.df)[-1], names_to = "Sample", values_to = "TPM") %>%
  mutate(Type = ifelse(grepl("-T-R" , Sample), "Tissue", "Organoid")) %>%
  mutate(log2TPM = log2(TPM+1)) %>%
  mutate(Gene_id = factor(Gene_id, levels = c("VEGFA", "FLT1", "KDR", "FLT4"))) %>%
  ggplot(aes(x=Type, y=log2TPM, fill=Type) ) + theme_base() +
  geom_dotplot(binaxis = "y", stackdir = "center", binpositions="all", binwidth = 0.3,color="black") +
  facet_grid(~ Gene_id)+ labs(y="log2 TPM", x =NULL) +
  theme(axis.text.x=element_text(color="black", angle = 45, hjust=1))
dev.off()



# Read TPM expression data and set gene IDs as row names
GBO_TPM <- read.delim("GBM_GBO_TPM.txt")
rownames(GBO_TPM) <- GBO_TPM$Gene_id

# Keep only columns with "OR" in their names (GBO samples)
GBO_TPM <- GBO_TPM[,grepl("OR", colnames(GBO_TPM))]

# Filter out genes with low average expression (<5 TPM)
GBO_TPM <- GBO_TPM[rowMeans(GBO_TPM) >= 5, ]

# Filter out genes where ≤3 samples have expression >1
GBO_TPM <- GBO_TPM %>%
  filter(rowSums(. > 1) > 3)

# Log2 transform the data (adding 1 to avoid log(0))
log_GBO_TPM.raw <- log2(GBO_TPM +1)
# Save current row names as a new column for filtering
log_GBO_TPM.raw$rownames <- rownames(log_GBO_TPM.raw)

# Keep genes with expression >0 in more than 3 samples
log_GBO_TPM <- log_GBO_TPM.raw %>%
  filter(rowSums(across(-rownames, ~ . > 0)) > 3)

# Restore original row names after filtering
rownames(log_GBO_TPM) <- rownames(GBO_TPM)

### Separate samples into TMZ-resistant and TMZ-sensitive groups

# TMZ-resistant samples
log_GBO_TPM.TMZ.resi <- log_GBO_TPM[,c("GBO.002.OR5.R", "GBO.003.OR4.R", "GBO.006.OR7.R")]
rownames(log_GBO_TPM.TMZ.resi) <- rownames(log_GBO_TPM)

# TMZ-sensitive samples
log_GBO_TPM.TMZ.sens <- log_GBO_TPM[,c("GBO.001.OR3.R", "GBO.005.OR4.R", "GBO.008.OR.p.3.R")]
rownames(log_GBO_TPM.TMZ.sens) <- rownames(log_GBO_TPM)


### Calculate fold-change between resistant and sensitive samples

log_GBO_TPM <- log_GBO_TPM %>% select(-rownames)

# Compute average expression per gene for each group
resi.mean <- rowMeans(log_GBO_TPM.TMZ.resi)
sens.mean <- rowMeans(log_GBO_TPM.TMZ.sens)
fold.diff <- resi.mean - sens.mean

# Select top 50 upregulated and top 50 downregulated genes (total 100)
top_DEGs <- c(names(sort(fold.diff, decreasing=TRUE)[1:50]), names(sort(fold.diff)[1:50]))


# Prepare data for the heatmap (excluding the last column if needed)
data <- log_GBO_TPM[top_DEGs, -ncol(log_GBO_TPM)]

# Rename columns for clarity
colnames(data) <- c("GBO-001", "GBO-002", "GBO-003", "GBO-005", "GBO-006","GBO-008")

# Define sample types for each column
sample_types <- c("TMZ-sensitive", "TMZ-resistant", "TMZ-resistant", 
                  "TMZ-sensitive", "TMZ-resistant", "TMZ-sensitive")

# Set colors for sample types
type_colors <- c("TMZ-sensitive" = "yellow", "TMZ-resistant" = "purple")

# Create column annotation data frame
annotation <- data.frame(Type = sample_types)
rownames(annotation) <- colnames(data)


# Define annotation colors for both column and row annotations
ann_colors <- list(
  Type = c("TMZ-sensitive" = "yellow", "TMZ-resistant" = "purple"), 
  `TMZ-resist Expr` = c("Upregulated" = "brown", "Downregulated" = "darkgreen") 
)

# Create row annotation indicating whether genes are up- or downregulated
annotation_row <- data.frame(`TMZ-resist Expr` = factor(c(rep("Upregulated", 50), rep("Downregulated", 50))))
rownames(annotation_row) <- rownames(data)
colnames(annotation_row) <- "TMZ-resist Expr"

# Plot the heatmap with row scaling and clustering
pheatmap(
  as.matrix(data),
  annotation_col = annotation,
  annotation_colors = ann_colors,
  annotation_row = annotation_row,
  annotation_names_row = FALSE,
  scale = "row", # Standardize by rows
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color=colorRampPalette(c("blue", "white", "red"))(100),
  show_rownames=TRUE
)



#Figure 4C, Figure S6

# Load TPM expression data and set row names to gene IDs
GBO_TPM <- read.delim("GBM_GBO_TPM.txt")
rownames(GBO_TPM) <- GBO_TPM$Gene_id
# Keep only columns with "OR" in their names (GBO samples)
GBO_TPM <- GBO_TPM[,grepl("OR", colnames(GBO_TPM))]

# Log2 transform the TPM data (adding 1 to avoid log(0))
log_GBO_TPM <- log2(GBO_TPM +1)

### Separate samples into TMZ-resistant and TMZ-sensitive groups
log_GBO_TPM.TMZ.resi <- log_GBO_TPM[,c("GBO.002.OR5.R", "GBO.003.OR4.R", "GBO.006.OR7.R")]
log_GBO_TPM.TMZ.sens <- log_GBO_TPM[,c("GBO.001.OR3.R", "GBO.005.OR4.R", "GBO.008.OR.p.3.R")]

# Calculate fold changes between resistant and sensitive samples
resi.mean <- rowMeans(log_GBO_TPM.TMZ.resi)
sens.mean <- rowMeans(log_GBO_TPM.TMZ.sens)
fold.diff <- resi.mean - sens.mean

# Print average expression for each group (for reference)
mean(resi.mean)
mean(sens.mean)

# Plot a histogram of the fold changes
hist(fold.diff, breaks=200, xlim=c(-3,3), xlab="Fold Change", main = "Histogram of GBO RNA expression FC")

# Calculate p-value for all genes (estimating differential expressions)
# Welch Two Sample t-test

# Calculate p-values using the Wilcoxon test for each gene
p.val <- sapply(1:nrow(log_GBO_TPM.TMZ.resi), function(i) 
  wilcox.test(as.matrix(log_GBO_TPM.TMZ.resi[i,]), 
              as.matrix(log_GBO_TPM.TMZ.sens[i,]), paired = FALSE)$p.value)

# Adjust p-values using Benjamini-Hochberg FDR
FDR <- round(p.adjust(p.val, 'BH'),3)

fData <- data.frame(Accession=rownames(log_GBO_TPM.TMZ.resi), FC=fold.diff, p.value=p.val, FDR=FDR)
fData

### Gene Ontology Enrichment Analysis

# Identify differentially expressed genes (DEGs)
fData$DEG.UP <- fData$FC > 1
fData$DEG.DOWN <- fData$FC < -1

# Select only DEGs (either up or down)
DEG <- fData[fData$DEG.UP | fData$DEG.DOWN,]
# Order DEGs by fold change (highest first)
DEG <- DEG[order(DEG$FC, decreasing=T),]

# Convert gene symbols to Entrez IDs using org.Hs.eg.db
entrezID <- bitr(DEG$Accession, fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

# Create a named vector of fold changes with Entrez IDs as names
FC <- DEG$FC
names(FC) <- entrezID$ENTREZID
FC

# Perform GO enrichment analysis (for all ontologies)
ego <- enrichGO(entrezID$ENTREZID,
                OrgDb = "org.Hs.eg.db",
                ont = "ALL",
                readable = TRUE,
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                universe = FC)

# Plot GO enrichment results with a dot plot, split by ontology
pdf("GO_dotplot_TMZresponse.pdf", height=16, width=12)
dotplot(ego, split = "ONTOLOGY", showCategory = 20, font.size = 12, label_format = 50) +
  facet_grid(ONTOLOGY ~ ., scale = "free_y") +
  geom_count() +
  scale_size_area(max_size = 8)
dev.off()


# Perform KEGG pathway enrichment analysis
kk <- enrichKEGG(entrezID$ENTREZID,
                 organism = 'hsa',
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05)

pdf("KEGG_dotplot_TMZresponse.pdf", height=16, width=12)
# Plot KEGG enrichment results using a dot plot
dotplot(kk, font.size = 10, showCategory = 20) +
  geom_count() +
  scale_size_area(max_size = 5)
dev.off()
