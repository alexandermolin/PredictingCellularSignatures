library(Seurat)
library(SeuratData)
library(ape)
library(dplyr)
library(Seurat)
library(patchwork)
# Load library for DESeq2
library(DESeq2)
# Load library for RColorBrewer
library(RColorBrewer)
# Load library for pheatmap
library(pheatmap)
# Load library for tidyverse
library(tidyverse)


# read colnames from txt:
colnames_from_txt <- function(name_of_file){
  return(scan(file = name_of_file, what = character(), nlines = 1, sep = "\t"))
}

# read rownames from txt:
rownames_from_txt <- function(name_of_file){
  return(scan(file = name_of_file, what = character(), flush = TRUE, sep = "\t"))
}

# read data from txt:
data_from_txt <- function(name_of_file, names_of_columns, names_of_rows){
  return(matrix(scan(file = name_of_file, skip = 0, what = character(), sep = "\t"),
                ncol = length(names_of_columns), nrow = length(names_of_rows),
                byrow = TRUE, dimnames = list(names_of_rows, names_of_columns)))
}

# read data from txt to sparse matrix with dimnames
txt_to_dgCMatrix <- function(name_of_file){
  spaltenbeschriftung <- colnames_from_txt(name_of_file)
  reihenbeschriftung <- rownames_from_txt(name_of_file)
  raw_data_read <- as.sparse(data_from_txt(name_of_file, spaltenbeschriftung, reihenbeschriftung))
  raw_data_read <- raw_data_read[, -1]
  raw_data_read <- raw_data_read[-1,]
  return(raw_data_read)
}

#SEurat####
## Load the PBMC dataset#####

# Initialize the Seurat object with the raw (non-normalized data).



filename_data_fractions_adjusted <- paste0("/home/eva/Documents/Alexander/endegelaende/test/cs_out/CIBERSORTxHiRes_NA_CTB_8W_2_Window13.txt")
stb_from_cs <- as.data.frame(txt_to_dgCMatrix(filename_data_fractions_adjusted))
stb_from_cs[is.na(stb_from_cs)] <- 0

#pbmc <- stb_from_cs

pbmc <- CreateSeuratObject(counts = stb_from_cs)

pbmc

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


##identification of highly variable features: ####
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

##pca####
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 10)

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")



DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)




# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:10)
JackStrawPlot(pbmc, dims = 1:10)

ElbowPlot(pbmc)

##cluster cells####
pbmc <- FindNeighbors(pbmc, dims = 1:7)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:5,  verbose = FALSE) #da kommt ein fehler..

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")


##find marker#####
# find all markers
if (requireNamespace("ape", quietly = TRUE)) {
  data(pbmc)
  pbmc
  pbmc <- BuildClusterTree(object = pbmc)
  Tool(object = pbmc, slot = 'BuildClusterTree')
}



pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)




#deSeq#####


filename_data_fractions_adjusted <- paste0("/home/eva/Downloads/CIBERSORTx_Job23_output/CIBERSORTxHiRes_Job23_CTB_Window13.txt")
stb_from_cs <- as.data.frame(txt_to_dgCMatrix(filename_data_fractions_adjusted))
stb_from_cs[is.na(stb_from_cs)] <- 0

sample_number_list_healthy <- c("SRR12363244", "SRR12363245", "SRR12363246", 
                                "SRR12363247", "SRR12363248", "SRR16760002",
                                "SRR16760003", "SRR16760004", "SRR16760005",
                                "SRR16760006", "SRR16760007", "SRR16760008",
                                "SRR16760009", "SRR16760010", "SRR16760012")

sample_number_list_unhealthy <- c("SRR13632931", "SRR13632932", "SRR13632933",
                                  "SRR13632934", "SRR13632935", "SRR14788891",
                                  "SRR14788893", "SRR14788895", "SRR14788897",
                                  "SRR14788899" )

# Create matrix


metadata_deseq <- as.data.frame(colnames(stb_from_cs))
rownames(metadata_deseq) <- metadata_deseq$`colnames(stb_from_cs)`
metadata_deseq$Health <- c(rep("healthy",length(metadata_deseq[,1])))
for (i in rownames(metadata_deseq)) {
  if (i %in% sample_number_list_unhealthy) {
    metadata_deseq$Health[rownames(metadata_deseq)==i] <- "unhealthy"
  }
}

metadata_deseq$Health <- as.factor(metadata_deseq$Health)


match(rownames(metadata_deseq), colnames(stb_from_cs))

idx <- match(rownames(metadata_deseq), colnames(stb_from_cs))

reordered_metaData <- metadata_deseq[idx,]

# check the order

all(rownames(reordered_metaData) == colnames(stb_from_cs))



dds <- DESeqDataSetFromMatrix(countData= round(stb_from_cs),
                              colData = metadata_deseq,
                              design = ~ 1)
dds <- estimateSizeFactors(dds)

sizeFactors(dds)
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))
normlzd_dds <- counts(dds, normalized=T)

head(normlzd_dds)
plot(hclust(dist(t(normlzd_dds))), labels=colData(dds)$protocol)
plot(hclust(dist(t(normlzd_dds))), labels=colData(dds)$Time)
plot(log(normlzd_dds[,1])+1, log(normlzd_dds[,2])+1, cex =.1)

# Varaiance Stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind = T)#, fitType = "mean")

# extract the vst matris from the object
vsd_mat <- assay(vsd)

# compute pairwise correlation values
vsd_cor <- cor(vsd_mat)

vsd_cor
pheatmap(vsd_cor)
plotPCA(vsd, intgroup = "Health")

# Calculating mean for each gene
mean_readCounts <- apply(stb_from_cs[,1:25], 1, mean)

# Calculating variance for each gene
var_readCounts <- apply(stb_from_cs[,1:25], 1, var)

df <- data.frame(mean_readCounts, var_readCounts)

# ggplot2 library
library(ggplot2)

ggplot(df) +
  geom_point(aes(x=mean_readCounts, y= var_readCounts)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Mean counts per gene") +
  ylab("Variance per gene") +
  labs(title = "DESeq2 model - Dispersion")


design(dds) <- ~ Health

dds <- DESeq(dds, fitType = "mean")
res <-results(dds)


summary(res)

plotMA(res, ylim=c(-5,5) )
resBigFC <- results(dds, lfcThreshold=1, altHypothesis="greaterAbs")
plotMA(resBigFC, ylim=c(-5,5))
abline(h=c(-1,1),lwd=5)
plotDispEsts(dds)
resSort <- res[order(res$pvalue),]

head(resSort)

##### ab da gehts nicht mehr, ##################

significantly_deg <- subset(resSort, resSort$pvalue < 0.05)
significantly_deg_2 <- data.frame(gene_name = rownames(significantly_deg)) 
significantly_deg_2$pVal <- significantly_deg$pvalue

res_sig <- data.frame(normlzd_dds[significantly_deg_2$gene_name, ])

table_for_print <- pivot_longer(significantly_deg, -c(gene_name),names_to = "unf", values_to = "PValues")

table_for_print <- table_for_print %>%
  mutate_if(is.numeric,
            round,
            digits = 5)


plot_table <- ggplot(table_for_print, aes(x=unf, y=gene_name, fill=PValues))+
  geom_tile(color = "black", lwd = 1.5, linetype = 1) + 
  geom_text(aes(label = PValues), color = "white", size = 4)+
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20, title = "P-Value"))+
  ylim(rev(levels(as.factor(table_for_print$gene_name)))) + xlim(levels(as.factor(table_for_print$unf)))

(plot_table + scale_fill_viridis_c(option = "D") + labs(x = "viridis", y = NULL))


heat_colors <- brewer.pal(6, "YlOrRd")

# Run pheatmap
pheatmap(res_sig,
         color = heat_colors,
         cluster_rows = T,
         show_rownames = T,
         annotation = dplyr::select(metadata_deseq, Health),
         scale = "row"
)


top_20 <- data.frame(normlzd_dds)[1:20, ] %>%
  rownames_to_column(var = "ensgene")

top_20 <- gather(top_20, 
                 key = "samplename", value = "normalized_counts", 2:8)

top_20 <- inner_join(top_20, rownames_to_column(reordered_metaData,
                                                var  = "samplename"), by = "samplename")

ggplot(top_20) +
  geom_point(aes(x = ensgene, y = normalized_counts, color = Health)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

dev.off()


save_mixture_file_for_cs(significantly_deg)

# 
# pdf(file = "fraction_plots.pdf",width = 14, height = 14)
# boxplot_celltypes
# density_binline_celltypes
# density_celltypes
# dev.off()

rm(list = ls())
gc()
