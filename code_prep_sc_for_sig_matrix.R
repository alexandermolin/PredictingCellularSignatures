# library(readxl)
# #library(DESeq2)
# library(SeuratObject)
# library(Seurat)
# library(dplyr)
# library(patchwork)
# library(ggplot2)
# library(umap)
# library(reticulate)
# library(EnsDb.Hsapiens.v79)
# library(tidyr)
# #library(biomaRt)
# #library(lubridate)
# library(Rtsne)
# library(tidyverse)
# 
# ####functions####
# 
# # to plot graphics into file:
# plot_to_pdf <- function(plot_titel){
#   return(pdf(paste0(plot_titel,"_", gsub('.{4}$', '', filename), ".pdf"), width = 34, height = 12))
# }
# 
# # read colnames from txt:
# colnames_from_txt <- function(name_of_file){
#   return(scan(file = name_of_file, what = character(), nlines = 1, sep = "\t"))
# }
# 
# # read rownames from txt:
# rownames_from_txt <- function(name_of_file){
#   return(scan(file = name_of_file, what = character(), flush = TRUE, sep = "\t"))
# }
# 
# # read data from txt:
# data_from_txt <- function(name_of_file, names_of_columns, names_of_rows){
#   return(matrix(scan(file = name_of_file, skip = 0, what = character(), sep = "\t"),
#                 ncol = length(names_of_columns), nrow = length(names_of_rows),
#                 byrow = TRUE, dimnames = list(names_of_rows, names_of_columns)))
# }
# 
# # read data and transform from txt:
# trans_data_from_txt <- function(name_of_file, names_of_columns, names_of_rows){
#   return(matrix(scan(file = name_of_file, skip = 0, what = character(), sep = "\t"),
#                 ncol = length(names_of_rows), nrow = length(names_of_columns),
#                 byrow = FALSE, dimnames = list(names_of_columns, names_of_rows)))
# }
# 
# # substring:
# substr_right <- function(a_vector, number_of_symbols){
#   sapply(a_vector, function(x)
#     substr(x, (nchar(x)-number_of_symbols+1), nchar(x))
#   )
# }
# 
# substr_left <- function(a_vector, number_of_symbols){
#   sapply(a_vector, function(x)
#     substr(x, 1, number_of_symbols)
#   )
# }
# 
# create_mixture_matrix <- function(counts_filenames, mapping, data_file_list){
#   for (i in 1:length(counts_filenames)) {
#     
#     spaltenbeschriftung <- colnames_from_txt(counts_filenames[i])
#     reihenbeschriftung <- rownames_from_txt(counts_filenames[i])
#     raw_data_read <- as.data.frame(data_from_txt(counts_filenames[i], spaltenbeschriftung, reihenbeschriftung))
#     raw_data_read <- raw_data_read[, -1]
#     raw_data_read <- raw_data_read[-1,]
#     raw_data_read$tx.tx_id <- rownames(raw_data_read)
#     raw_data_read <- subset(raw_data_read, select= c(tx.tx_id, tpm))
#     
#     # Rownames to ensembl-ID
#     raw_data_read$tx.tx_id <- sapply(raw_data_read$tx.tx_id,function(x)
#       substr(x, 1, 15))
#     
#     mapping <- full_join(mapping, raw_data_read)
#     mapping$tpm <- replace(mapping$tpm, which(is.na(mapping$tpm)), 0)
#     colnames(mapping)[which(names(mapping) == "tpm")] <- substr_right(data_file_list[i], 11)
#   }
#   mapping$tx.tx_id <- NULL
#   mapping$tx.gene_id <- NULL
#   reihennamen_fuer_matrix <- mapping$tx.gene_name
#   mapping$tx.gene_name <- NULL
#   
#   for (i in colnames(mapping[1:length(colnames(mapping))])) {
#     mapping[,i] <- as.numeric(mapping[,i])
#   }
#   
#   mapping_matrix <- as.matrix(mapping)
#   rownames(mapping_matrix) <- reihennamen_fuer_matrix
#   
#   mapping_matrix_neu <- t(sapply(by(mapping_matrix,rownames(mapping_matrix),colSums),identity))
#   return(mapping_matrix_neu)
# }
# 
# ### for manual run with timestamp:
# # save_mixture_file_for_cs <- function(mixture_matrix){
# #   write.table(as.matrix(mixture_matrix), sep = "\t", file = gsub( ":", "_", paste0(substitute(mixture_matrix), as.ITime(Sys.time()), ".txt" )), quote = FALSE,
# #               row.names = TRUE, col.names = NA)
# # }
# 
# save_mixture_file_for_cs <- function(mixture_matrix){
#   write.table(as.matrix(mixture_matrix), sep = "\t", file = paste0(substitute(mixture_matrix), ".txt" ), quote = FALSE,
#               row.names = TRUE, col.names = NA)
# }
# 
# save_gene_subset_file_for_cs <- function(gene_subset_file){
#   write.table(as.matrix(gene_subset_file), sep = "\t", file = paste0(substitute(gene_subset_file), ".txt" ), quote = FALSE,
#               row.names = FALSE, col.names = FALSE)
# }
# 
# 
# 
# get_script_path <- function(){
#   cmd.args <- commandArgs()
#   m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
#   script_dir <- dirname(regmatches(cmd.args, m))
#   if(length(script_dir) == 0) stop("can't determine script dir: please call the script with Rscript")
#   if(length(script_dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
#   return(script_dir)
# }
# 
# # read data from txt to sparse matrix with dimnames
# txt_to_dgCMatrix <- function(name_of_file){
#   spaltenbeschriftung <- colnames_from_txt(name_of_file)
#   reihenbeschriftung <- rownames_from_txt(name_of_file)
#   raw_data_read <- as.sparse(data_from_txt(name_of_file, spaltenbeschriftung, reihenbeschriftung))
#   raw_data_read <- raw_data_read[, -1]
#   raw_data_read <- raw_data_read[-1,]
#   return(raw_data_read)
# }
# 
# 
# 
# ####directory,files, etc####
# 
# 
# 
# #setting workingdirectory 
# ### for manual run:
# setwd("~/Documents/Alexander/endegelaende/test/")
# 
# #setwd(get_script_path())
# 
# 
# 
# 
# ####sc_signature_matrix:#####
# 
# 
# #data read:
# data_file_name <- "/home/eva/Documents/Alexander/endegelaende/test/c_single_cell_rna_data_hongmei.txt"
# raw_data_read <- txt_to_dgCMatrix(data_file_name) 
# 
# data_file_name <- "/home/eva/Documents/Alexander/endegelaende/test/c_tsne_df_von_xls.csv"
# spaltenbeschriftung <- colnames_from_txt(data_file_name)
# reihenbeschriftung <- rownames_from_txt(data_file_name)
# c_tsne_df <- as.data.frame(data_from_txt(data_file_name, spaltenbeschriftung, reihenbeschriftung))
# c_tsne_df <- c_tsne_df[, -1]
# c_tsne_df <- c_tsne_df[-1,]
# 
# 
# raw_data_read_orig <- raw_data_read
# 
# for (j in 1:ncol(raw_data_read)) {
#   colnames(raw_data_read)[j] <- c_tsne_df[colnames(raw_data_read)[j], "Apcluster"]
# }
# 
# raw_data_read <- raw_data_read[,(colnames(raw_data_read) %in% c_tsne_df$Apcluster)]
# #raw_data_read[1:4,1:4]
# raw_data_missing <- raw_data_read[,!(colnames(raw_data_read) %in% c_tsne_df$Apcluster)]
# #as.factor(raw_data_read)
# 
# #check geneexpression !=0
# 
# #raw_data_read <- raw_data_read[,colnames(raw_data_read) [colSums(raw_data_read)!=0]]
# any(!(colSums(raw_data_read) > 0))
# sc_for_signature_matrix <- raw_data_read[,!(colnames(raw_data_read) %in% c("EVT_24W_1", "EVT_24W_2"))]
# 
# raw_data_read[1:10,1:18]
# raw_data_read_orig[1:10,1:18]
# sc_for_signature_matrix[1:10,1:18]
# 
# ####sc_signature_matrix_merged:#####
# ####beschriftungen:
# 
# merged_rules <- unique(data.frame(genes = c_tsne_df$Apcluster, class = c_tsne_df$Merge_Type))
# 
# 
# merged_classes <- as.data.frame(colnames(raw_data_read))  
# colnames(merged_classes) <- "genes"
# 
# merged_classes <- left_join(merged_classes, merged_rules)
# sum(summary(as.factor(merged_classes$class)))
# merged_classes$genes <- NULL
# 
# save_gene_subset_file_for_cs(merged_classes)
# 
# sc_for_signature_matrix_merged <- raw_data_read
# colnames(sc_for_signature_matrix_merged) <- merged_classes$class
# 
# ####nonsense matrix:####
# werte_fuer_nonsense <- as.vector(raw_data_read)
# set.seed(1)
# sc_for_signature_matrix_merged_max_nonsense <- as.sparse(matrix(sample(werte_fuer_nonsense, length(werte_fuer_nonsense), replace = FALSE),
#                                                       nrow = 20866, ncol = 1449, dimnames = list(
#                                                         rownames(sc_for_signature_matrix_merged),
#                                                         colnames(sc_for_signature_matrix_merged))))
# nonsensecolnames <- rep_len(LETTERS[1:4], length(colnames(sc_for_signature_matrix_merged_max_nonsense)))
# colnames(sc_for_signature_matrix_merged_max_nonsense) <- nonsensecolnames
# 
# rm(merged_rules, merged_classes, raw_data_missing, nonsensecolnames, reihenbeschriftung, 
#    spaltenbeschriftung, werte_fuer_nonsense, data_file_list_unhealthy, data_file_list_healthy, data_file_list,
#    data_file_name)
# 
# #####seurat####
# #fuer_pbmc <- raw_data_read_orig
# #fuer_pbmc <- raw_data_read
# fuer_pbmc <- sc_for_signature_matrix
# #fuer_pbmc <- sc_for_signature_matrix_merged
# #fuer_pbmc <- sc_for_signature_matrix_merged_max_nonsense
# 
# #colnames(fuer_pbmc) <- gsub("_", "", colnames(fuer_pbmc))
# 
# cbmc <- CreateSeuratObject(counts = fuer_pbmc,
#                           # min.cells = 3, min.features = 200,
#                            names.field = 1,
#                           # names.delim = ".",
#                           # meta.data = c_tsne_df,
#                           row.names = rownames(fuer_pbmc))
# cbmc
# 
# set.seed(42)
# # perform visualization and clustering steps
# cbmc <- NormalizeData(cbmc)
# cbmc <- FindVariableFeatures(cbmc)
# cbmc <- ScaleData(cbmc)
# cbmc <- RunPCA(cbmc, verbose = TRUE)
# #ElbowPlot(cbmc)
# cbmc <- FindNeighbors(cbmc, dims = 1:50)
# cbmc <- FindClusters(cbmc, resolution = 0.8, verbose = TRUE) #resolution: 0.025 fuer 7 cluster, 1.09 fuer 14 cluster, 0.724 fuer 1 bei nonsense
# cbmc <- RunUMAP(cbmc, dims = 1:50)
# cbmc <- RunTSNE(cbmc, dims = 1:50)
# 
# 
# 
# #umap:
# cbmc_plot <- DimPlot(cbmc, label = TRUE)+
#   ggtitle("UMAP plot - sig.matrix")
# cbmc_plot
# 
# # png(filename = paste0("umap_sc_","sig_mat",".png"))
# # cbmc_plot
# # dev.off()
# 
# # #tsne:
# # cbmc[["umap"]] <- NULL
# # DimPlot(cbmc, label = TRUE)
# 
# # #PCA:
# # cbmc[["tsne"]] <- NULL
# # DimPlot(cbmc, label = TRUE)
# # 
# 
# # 
# # png(filename = paste0("cbmc_sc","sig_nonsense",".png"))
# # cbmc_plot
# # dev.off()
# 
# # cbmc.no.umap <- cbmc
# # cbmc.no.umap[["tsne"]] <- NULL
# # DimPlot(cbmc.no.umap) + RotatedAxis()
# 
# 
# # 
# # ###tsne:
# # seurat <- FindClusters(cbmc, resolution = seq(from = 0.1, to = 1, by = 0.1))
# # ResolutionList <- grep("_snn_res", colnames(seurat@meta.data), value = TRUE)
# # 
# # for (Resolution in ResolutionList){
# #   pdf(paste0(Resolution, "_tSNE.pdf"), width=7, height=7)
# #   g <- TSNEPlot(object = seurat, group.by = Resolution)
# #   print(g)
# #   dev.off()
# # }
# # 
# # Dim
# # pff_matrix <- as.data.frame(sc_for_signature_matrix)
# # # 
# # # umap_fit_TP <- umap(pff_matrix)
# # # 
# # #UMAP dataframe creation
# # umap_df_TP <- as.data.frame(umap_fit_TP$layout) 
# # names(umap_df_TP) <- c("UMAP1", "UMAP2")
# # umap_df_TP <- cbind(umap_df_TP, meta.TP.full)
# # umap_df_TP$Replicate <- as.factor(umap_df_TP$Replicate)
# # 
# # #UMAP plot
# # ggplot(umap_df_TP, aes(x=UMAP1,y=UMAP2, label="Name")) +
# #   scale_colour_brewer(palette="Dark2") +
# #   geom_point(size=4) + geom_text(vjust=2) +
# #   theme_bw(base_size = 22) + 
# #   theme(legend.position = "right") + #without legend type "none"
# #   ggtitle("UMAP plot - TP - all metabolites")
# 
# 
# 
# # ###tsne mit den coordinaten aus dem Paper metafile:
# # c_tsne_df$tSNE1 <- as.numeric(c_tsne_df$tSNE1)
# # c_tsne_df$tSNE2 <- as.numeric(c_tsne_df$tSNE2)
# # ggplot(c_tsne_df, aes(x = tSNE1, y = tSNE2, label = "")) +
# #   geom_point(size = 1) + geom_text(vjust=1) +
# #   theme_bw(base_size = 22) + 
# #   theme(legend.position = "none") + #without legend type "none"
# #   ggtitle("tSNE vom Paper")
# # 
# 
# # 
# # e combined several methods for cell sample clustering. The
# # dissimilarity matrix was built by calculating (1-β)/2, where β is the
# # Spearman correlation of each cell pair. Hierarchical clustering was
# # performed on the dissimilarity matrix using the ‘hclust’ R function
# # with the ‘average’ method, and the clusters were identified using
# # the R function ‘cutreeDynamic’ with the ‘hybrid’ method.58 The
# # clusters and cell types were visualized via t-distributed stochastic
# # neighbor embedding (t-SNE, as implemented in the’Rtsne’
# #                     package) on the dissimilarity matrix. PCA was performed using
# # the Seurat package with highly variant genes.
# 
# # hclust(pff_matrix, method = "average", members = NULL)
# 
# 
# 
# ####speichern####
# 
# save_mixture_file_for_cs(sc_for_signature_matrix_merged_max_nonsense)
# save_mixture_file_for_cs(sc_for_signature_matrix_merged)
# save_mixture_file_for_cs(sc_for_signature_matrix)
# 
# sc_for_signature_matrix_merged_max_nonsense[1:8,1:32]
# sc_for_signature_matrix_merged[1:8,1:32]
# raw_data_read[1:8,1:32]
# 
# rm(list = ls())
# gc()
# 
# 
