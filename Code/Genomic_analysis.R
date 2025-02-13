# genomic analysis

library(maftools)
library(tidyverse)
library(VennDiagram)
library(reshape2)
library(gridExtra)
library(ggpubr)
library(ggthemes)
library(ggExtra)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)


# oncoplot
anno = read.delim("anno.txt", row.names = 1)

gbo_paired_maf <- read.delim("GBO_paired.maf.txt", sep="\t")

gbo_paired_maf$Patient <- substr(gbo_paired_maf$Tumor_Sample_Barcode, 1 ,7) %>% gsub("GBO-","P-",.) %>% gsub("P-019","P-005",.)
gbo_paired_maf$Sample <- substr(gbo_paired_maf$Tumor_Sample_Barcode, 1 ,7)
gbo_paired_maf$Type <- grepl("T-D",gbo_paired_maf$Tumor_Sample_Barcode) %>% gsub("TRUE","Tissue",.) %>% gsub("FALSE","Organoid",.)

gbo_paired_maf_final <- gbo_paired_maf %>% filter(VAF >=0.05)
paired_maf <- read.maf(gbo_paired_maf_final, clinicalData=anno)

maf_col = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(maf_col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 
                   'Multi_Hit', 'Frame_Shift_Ins',
                   'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')

pdf("oncoplot.pdf", 7, 10)
oncoplot( maf = paired_maf , showTumorSampleBarcodes = T, 
          clinicalFeatures=c('Patient', 'Type', 'Extent_of_Resection', 'MGMT_methylation'), 
          sampleOrder = c("GBO-002-T-D", "GBO-002-OR5-D","GBO-003-T-D" ,"GBO-003-OR4-D" ,"GBO-005-T-D", "GBO-005-OR4-D", "GBO-005R-OR2-D" ,"GBO-008-T-D","GBO-008-OR-p-3-D"  ), 
          colors=maf_col, barcode_mar = 8, removeNonMutated = F, 
          draw_titv=T, gene_mar = 7 , top = 30,
          annotationColor = list(Extent_of_Resection = c("GTR" = "#F24738", "NTR" = "#F29B30", "PR" = "#328BD9", "STR" = "#025E73"), 
                                 MGMT_methylation = c("1" = "#9498F2", "0" = "#84BF04"), 
                                 Patient = c("P-002"="#BF634E","P-003"="#5D568C","P-005"="#F2F0CE","P-008"="#F27127"),
                                 Type = c("Tissue" = "#B3CDE3", "Organoid" = "#FBB4AE")),
          annoBorderCol = "white", sepwd_samples = 1, showTitle = F, anno_height = 1.5, SampleNamefontSize = 1.2)
dev.off()

# venn diagram
gbo_paired_maf_final$Mut_full <- paste(gbo_paired_maf_final$Hugo_Symbol,gbo_paired_maf_final$Chromosome,gbo_paired_maf_final$Start_Position,gbo_paired_maf_final$Variant_Classification,gbo_paired_maf_final$HGVSp_Short,sep="_")

ven_table<-gbo_paired_maf_final %>% filter(Patient == "P-003") %>% select("Mut_full","Tumor_Sample_Barcode") %>% table %>% as.data.frame.matrix() %>% table %>% as.data.frame.matrix()
num_area1<-ven_table[2,1]+ven_table[2,2] # organoid sum
num_area2<-ven_table[1,2]+ven_table[2,2] # tissue sum
num_cross<-ven_table[2,2]
plot<-draw.pairwise.venn(area1 = num_area1, area2 = num_area2, cross.area = num_cross,
                           fill=c("#FBB4AE","#B3CDE3"),cex = rep(6, 3),lwd=c(3,3),ext.text=TRUE,margin=0.1,
                           inverted= (num_area2>num_area1),ext.pos=90,ext.length=c(0),ext.dist=c(0,0))
assign("P003_ven",plot)

ven_table<-gbo_paired_maf_final %>% filter(Sample == "GBO-005") %>% select("Mut_full","Tumor_Sample_Barcode") %>% table %>% as.data.frame.matrix() %>% table %>% as.data.frame.matrix()
num_area1<-ven_table[2,1]+ven_table[2,2] # organoid sum
num_area2<-ven_table[1,2]+ven_table[2,2] # tissue sum
num_cross<-ven_table[2,2]
plot<-draw.pairwise.venn(area1 = num_area1, area2 = num_area2, cross.area = num_cross,
                           fill=c("#FBB4AE","#B3CDE3"),cex = rep(6, 3),lwd=c(3,3),ext.text=TRUE,margin=0.1,
                           inverted= (num_area2>num_area1),ext.pos=90,ext.length=c(0),ext.dist=c(0,0))
assign("P005_ven",plot)

example_OT.ven<-draw.pairwise.venn(area1 = 50, area2 = 50, cross.area = 25,
                                   fill=c("#FBB4AE","#B3CDE3"),cat.cex = rep(5, 2),lwd=c(3,3),margin=0.1,cex=c(0),
                                   ext.pos=90,ext.length=c(0),ext.dist=c(0,0),category = c("Organoid","Tissue"),cat.dist=0.1,cat.pos=c(-30,30))

ven_table<-gbo_paired_maf_final %>% filter(Patient == "P-008") %>% select("Mut_full","Tumor_Sample_Barcode") %>% table %>% as.data.frame.matrix() %>% table %>% as.data.frame.matrix()
num_area1<-ven_table[1,1]+ven_table[1,2] # organoid sum
num_area2<-ven_table[1,2] # tissue sum
num_cross<-ven_table[1,2]
plot<-draw.pairwise.venn(area1 = num_area1, area2 = num_area2, cross.area = num_cross,
                           fill=c("#FBB4AE","#B3CDE3"),cex = rep(6, 3),lwd=c(3,3),ext.text=TRUE,margin=0.1,
                           inverted= (num_area2>num_area1),ext.pos=90,ext.length=c(0),ext.dist=c(0,0))
assign("P008_ven",plot)

pdf("Venn.pdf",8,8)
grid.arrange(gTree(children=P003_ven), top=textGrob("P-003", gp=gpar(fontsize=50)))
grid.arrange(gTree(children=P005_ven), top=textGrob("P-005", gp=gpar(fontsize=50)))
grid.arrange(gTree(children=P008_ven), top=textGrob("P-008", gp=gpar(fontsize=50)))
grid.arrange(gTree(children=example_OT.ven))
dev.off()

# VAF correlation
gbo_paired_maf_final %>% 
  filter(Sample %in% c("GBO-003","GBO-005","GBO-008")) %>% 
  select("Sample", "Tumor_Sample_Barcode", "Type", "Mut_full", "VAF") %>% 
  dcast(Sample+Mut_full ~ Type, value.var="VAF") %>%
  drop_na(Tissue) %>% drop_na(Organoid) %>%
  ggplot(aes(x=Tissue, y=Organoid)) + geom_point(aes(col=Sample)) + stat_cor(method="spearman")+
  theme_base() + theme(axis.text.x=element_text(color="black",size=13),
                       axis.text.y=element_text(color="black",size=13),
                       axis.title=element_text(color="black", size=13)) +
  labs(x= "Tissue VAF", y="Organoid VAF") -> p
ggMarginal(p, type = "densigram", xparams = list(fill = "grey"),yparams = list(fill = "grey")) -> p
p
pdf("VAF_correlation.pdf", 6, 5)
p;dev.off()

# mutational siganture
#library(devtools)
#install_github("raerose01/deconstructSigs")
library(deconstructSigs)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)

gbo_paired_maf_final.sigdat <- gbo_paired_maf_final %>% filter(Variant_Type == "SNP") %>%
  dplyr::select("Tumor_Sample_Barcode","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2")

sigs.input <- mut.to.sigs.input(mut.ref = gbo_paired_maf_final.sigdat, 
                                bsg = BSgenome.Hsapiens.UCSC.hg38 ,
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_Position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2")

sig_GBO_002_OR5_D <- whichSignatures(tumor.ref = sigs.input,
  signatures.ref = signatures.exome.cosmic.v3.may2019,
  sample.id = 'GBO-002-OR5-D',contexts.needed = TRUE)
sig_GBO_003_OR4_D <- whichSignatures(tumor.ref = sigs.input,
  signatures.ref = signatures.exome.cosmic.v3.may2019,
  sample.id = 'GBO-003-OR4-D',contexts.needed = TRUE)
sig_GBO_005_OR4_D <- whichSignatures(tumor.ref = sigs.input,
  signatures.ref = signatures.exome.cosmic.v3.may2019,
  sample.id = 'GBO-005-OR4-D',contexts.needed = TRUE)
sig_GBO_008_OR_p_3_D <- whichSignatures(tumor.ref = sigs.input,
  signatures.ref = signatures.exome.cosmic.v3.may2019,
  sample.id = 'GBO-008-OR-p-3-D',contexts.needed = TRUE)

sig_GBO_002_T_D <- whichSignatures(tumor.ref = sigs.input,
  signatures.ref = signatures.exome.cosmic.v3.may2019,
  sample.id = 'GBO-002-T-D',contexts.needed = TRUE)
sig_GBO_003_T_D <- whichSignatures(tumor.ref = sigs.input,
  signatures.ref = signatures.exome.cosmic.v3.may2019,
  sample.id = 'GBO-003-T-D',contexts.needed = TRUE)
sig_GBO_005_T_D <- whichSignatures(tumor.ref = sigs.input,
  signatures.ref = signatures.exome.cosmic.v3.may2019,
  sample.id = 'GBO-005-T-D',contexts.needed = TRUE)
sig_GBO_008_T_D <- whichSignatures(tumor.ref = sigs.input,
  signatures.ref = signatures.exome.cosmic.v3.may2019,
  sample.id = 'GBO-008-T-D',contexts.needed = TRUE)

plotSignatures(sig_GBO_002_OR5_D)
makePie(sig_GBO_002_OR5_D)

plotSignatures(sig_GBO_003_OR4_D)
makePie(sig_GBO_003_OR4_D)

plotSignatures(sig_GBO_005_OR4_D)
makePie(sig_GBO_005_OR4_D)

plotSignatures(sig_GBO_008_OR_p_3_D)
makePie(sig_GBO_008_OR_p_3_D)

plotSignatures(sig_GBO_002_T_D)
makePie(sig_GBO_002_T_D)

plotSignatures(sig_GBO_003_T_D)
makePie(sig_GBO_003_T_D)

plotSignatures(sig_GBO_005_T_D)
makePie(sig_GBO_005_T_D)

plotSignatures(sig_GBO_008_T_D)
makePie(sig_GBO_008_T_D)


# VAF in Patient 5
gbo_paired_maf_final %>% filter(Tumor_Sample_Barcode %in% c("GBO-005-OR4-D", "GBO-019-OR2-D")) %>% 
  select(Tumor_Sample_Barcode, Mut_full) %>% 
  dcast(Mut_full ~ ., value.var="Tumor_Sample_Barcode",  fun.aggregate=function(x) paste(x, collapse = ", ")) -> PDO5_info

colnames(PDO5_info)[2] <- "Mut_in"

PDO5_maf <- gbo_paired_maf_final %>% filter(Tumor_Sample_Barcode %in% c("GBO-005-OR4-D", "GBO-019-OR2-D"))
PDO5_maf <- merge(PDO5_maf, PDO5_info, by="Mut_full")
PDO5_maf$Mut_in_tissue <- PDO5_maf$Mut_full %in% subset(gbo_paired_maf_final, Tumor_Sample_Barcode == "GBO-005-T-D")$Mut_full

pdf("GBO-005R_VAF_boxplot.pdf", 6,6)
PDO5_maf %>% filter(Tumor_Sample_Barcode == "GBO-019-OR2-D") %>% 
  filter(c(Mut_in == "GBO-005-OR4-D, GBO-019-OR2-D" & Mut_in_tissue == "TRUE") |
         c(Mut_in == "GBO-019-OR2-D" & Mut_in_tissue == "FALSE")) %>%
  mutate(Mut_in = gsub("GBO-005-OR4-D, GBO-019-OR2-D", "GBM-005 & GBO-005", Mut_in)) %>%
  mutate(Mut_in = gsub("GBO-019-OR2-D", "GBO-005R", Mut_in)) %>%
  mutate(Primary_tissue = ifelse(Mut_in_tissue, "Mutation", "No mutation")) %>%
  ggplot(aes(x=Mut_in, y=VAF)) + geom_boxplot(outlier.shape = NA, aes(fill=Mut_in)) + 
  geom_jitter(width=0.2, height = 0, size = 2, alpha = 0.8 , aes(shape=Primary_tissue)) +
  theme_bw() + theme(axis.text.x=element_text(color="black",size=13),
                     axis.text.y=element_text(color="black",size=13),
                     axis.title=element_text(color="black", size=13)) + 
  stat_compare_means(label.x = 1.8) + labs(x=NULL) +
  guides(fill = "none") +
  ggtitle("VAF in GBO-005R")
dev.off()


### Figure S9

# Load CNV data and filter for selected genes
GBO_CNV <- read.delim("GBM_GBO_CNV.txt")
Genes <- c(TME_related, PI3K_pathway, RAS_pathway)
GBO_CNV.f <- GBO_CNV[GBO_CNV$Gene %in% Genes,]

# Label samples as GBO or GBM based on the Sample name
GBO_CNV.f <- GBO_CNV.f %>%
  mutate(Type = if_else(grepl("OR", Sample), "GBO", "GBM"))

# Create a mechanism of action (MoA) annotation for each gene
MoA <- c(rep("TME maintenance", length(TME_related)), rep("PI3K/AKT/MTOR", length(PI3K_pathway)), rep("RAS/RAF/ERK", length(RAS_pathway)))
Genes <- c(TME_related, PI3K_pathway, RAS_pathway)
MoA_df <- data.frame(Genes, MoA)

# Mark samples with "normal" CNV as NA
GBO_CNV.f$Type[GBO_CNV.f$CNV == "normal"] <- "NA"


# Select columns and reshape data to wide format
long_data <- GBO_CNV.f[, c("Sample", "Gene", "Log2_Ratio")]
wide_data <- dcast(long_data, Gene ~ Sample, value.var = "Log2_Ratio", fun.aggregate = mean)
rownames(wide_data) <- wide_data$Gene
wide_data <- wide_data[, -1]  

# Rename columns for clarity
colnames(wide_data) <- c("GBO-001", "GBM-001", "GBO-002", "GBM-002", "GBO-003", "GBM-003", "GBO-005", "GBM-005", "GBO-006", "GBO-008", "GBM-008", "GBO-019")

# Define gene sets for each pathway(from MsigDB)
TME_related <- c("PDGFRB", "FGFR1", "FGFR2", "FGFR3", "FGFR4")
PI3K_pathway <- c("EGFR", "AKT1", "AKT2", "AKT3", "PIK3CA", "PIK3CB", "PIK3CD", "MTOR")
RAS_pathway <- c("GRB2", "SOS1", "SOS2", "KRAS", "HRAS", "NRAS", "RAF1", "ARAF", "BRAF", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2")

# Annotate each gene by its pathway
pathway_annotation <- ifelse(rownames(wide_data) %in% TME_related, "TME_related",
                             ifelse(rownames(wide_data) %in% PI3K_pathway, "PI3K_pathway", "RAS_pathway"))

# Choose colors for pathway annotation
mycol=brewer.pal(4,"Dark2")


# Create row annotation for the heatmap
row_anno <- rowAnnotation(Pathway = pathway_annotation,
                          col = list(Pathway = c("TME_related" = mycol[1], 
                                                 "PI3K_pathway" = mycol[2], 
                                                 "RAS_pathway" = mycol[3])))


# Determine the minimum and maximum values in the data (ignoring NA)
data_min <- min(wide_data, na.rm = TRUE)
data_max <- max(wide_data, na.rm = TRUE)


pdf("Heatmap_log2CNV.pdf", width=8, height=6)

Heatmap(as.matrix(wide_data), 
        name = "Log2 Ratio", 
        show_row_names = TRUE, 
        column_order = colnames(wide_data),
        row_order = c(PI3K_pathway, RAS_pathway, TME_related),
        show_column_names = TRUE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        right_annotation = row_anno,
        col = colorRamp2(c(data_min, 0, data_max), c("blue", "white", "red")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, 
                    gp = gpar(col = "grey", lwd = 0.5, fill = NA))
        })

dev.off()
