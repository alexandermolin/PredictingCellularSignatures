library(readxl)
#library(DESeq2)
library(SeuratObject)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(umap)
library(reticulate)
library(EnsDb.Hsapiens.v79)
library(tidyr)
#library(biomaRt)
#library(lubridate)
library(Rtsne)
library(tidyverse)

####functions####

# to plot graphics into file:
plot_to_pdf <- function(plot_titel){
  return(pdf(paste0(plot_titel,"_", gsub('.{4}$', '', filename), ".pdf"), width = 34, height = 12))
}

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

# read data and transform from txt:
trans_data_from_txt <- function(name_of_file, names_of_columns, names_of_rows){
  return(matrix(scan(file = name_of_file, skip = 0, what = character(), sep = "\t"),
                ncol = length(names_of_rows), nrow = length(names_of_columns),
                byrow = FALSE, dimnames = list(names_of_columns, names_of_rows)))
}

# substring:
substr_right <- function(a_vector, number_of_symbols){
  sapply(a_vector, function(x)
    substr(x, (nchar(x)-number_of_symbols+1), nchar(x))
  )
}

substr_left <- function(a_vector, number_of_symbols){
  sapply(a_vector, function(x)
    substr(x, 1, number_of_symbols)
  )
}

create_mixture_matrix <- function(counts_filenames, mapping, data_file_list){
  for (i in 1:length(counts_filenames)) {
    
    spaltenbeschriftung <- colnames_from_txt(counts_filenames[i])
    reihenbeschriftung <- rownames_from_txt(counts_filenames[i])
    raw_data_read <- as.data.frame(data_from_txt(counts_filenames[i], spaltenbeschriftung, reihenbeschriftung))
    raw_data_read <- raw_data_read[, -1]
    raw_data_read <- raw_data_read[-1,]
    raw_data_read$tx.tx_id <- rownames(raw_data_read)
    raw_data_read <- subset(raw_data_read, select= c(tx.tx_id, tpm))
    
    # Rownames to ensembl-ID
    raw_data_read$tx.tx_id <- sapply(raw_data_read$tx.tx_id,function(x)
      substr(x, 1, 15))
    
    mapping <- full_join(mapping, raw_data_read)
    mapping$tpm <- replace(mapping$tpm, which(is.na(mapping$tpm)), 0)
    colnames(mapping)[which(names(mapping) == "tpm")] <- substr_right(data_file_list[i], 11)
  }
  mapping$tx.tx_id <- NULL
  mapping$tx.gene_id <- NULL
  reihennamen_fuer_matrix <- mapping$tx.gene_name
  mapping$tx.gene_name <- NULL
  
  for (i in colnames(mapping[1:length(colnames(mapping))])) {
    mapping[,i] <- as.numeric(mapping[,i])
  }
  
  mapping_matrix <- as.matrix(mapping)
  rownames(mapping_matrix) <- reihennamen_fuer_matrix
  
  mapping_matrix_neu <- t(sapply(by(mapping_matrix,rownames(mapping_matrix),colSums),identity))
  return(mapping_matrix_neu)
}

### for manual run with timestamp:
# save_mixture_file_for_cs <- function(mixture_matrix){
#   write.table(as.matrix(mixture_matrix), sep = "\t", file = gsub( ":", "_", paste0(substitute(mixture_matrix), as.ITime(Sys.time()), ".txt" )), quote = FALSE,
#               row.names = TRUE, col.names = NA)
# }

save_mixture_file_for_cs <- function(mixture_matrix){
  write.table(as.matrix(mixture_matrix), sep = "\t", file = paste0(substitute(mixture_matrix), ".txt" ), quote = FALSE,
              row.names = TRUE, col.names = NA)
}

save_gene_subset_file_for_cs <- function(gene_subset_file){
  write.table(as.matrix(gene_subset_file), sep = "\t", file = paste0(substitute(gene_subset_file), ".txt" ), quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}



get_script_path <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script_dir <- dirname(regmatches(cmd.args, m))
  if(length(script_dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script_dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script_dir)
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



####directory,files, etc####



#setting workingdirectory 
### for manual run:
setwd("~/Documents/Alexander/endegelaende/test/")

#setwd(get_script_path())


####MAPPING BULK:
#Data-directory
#dirname <- "~/Documents/Alexander/endegelaende/data_bulk/"

#list of fastq-files
data_file_list <- list.files(all.files = FALSE, full.names = FALSE, pattern = "^output_kallisto_")


#all samples:
counts_filenames <- sapply(data_file_list, function(x)
                          paste0(x,"/abundance.tsv"))



#healthy samples:
data_file_list_healthy <- c("output_kallisto_SRR12363244", "output_kallisto_SRR12363245", "output_kallisto_SRR12363246",
                            "output_kallisto_SRR12363247", "output_kallisto_SRR12363248", "output_kallisto_SRR16760002",
                            "output_kallisto_SRR16760003", "output_kallisto_SRR16760004", "output_kallisto_SRR16760005",
                            "output_kallisto_SRR16760006", "output_kallisto_SRR16760007", "output_kallisto_SRR16760008",
                            "output_kallisto_SRR16760009", "output_kallisto_SRR16760010", "output_kallisto_SRR16760012")

counts_filenames_healthy <- sapply(data_file_list_healthy, function(x)
  paste0(x,"/abundance.tsv"))



data_file_list_unhealthy <- c("output_kallisto_SRR13632931", "output_kallisto_SRR13632932", "output_kallisto_SRR13632933",
                              "output_kallisto_SRR13632934", "output_kallisto_SRR13632935", "output_kallisto_SRR14788891",
                              "output_kallisto_SRR14788893", "output_kallisto_SRR14788895", "output_kallisto_SRR14788897",
                              "output_kallisto_SRR14788899" )

counts_filenames_unhealthy <- sapply(data_file_list_unhealthy, function(x)
  paste0(x,"/abundance.tsv"))


####create a df with all Genes####
edb <- EnsDb.Hsapiens.v79
## Get all transcripts defined in Ensembl (version 79):
tx <- transcripts(edb, columns=c("tx_id", "gene_id", "gene_name"))
## extract the transcript ids and gene names
mapping <- data.frame(tx$tx_id, tx$gene_name)
rownames(mapping) <- mapping$tx.tx_id
#dim(mapping)
#head(mapping)

#### Fill df with bulk data ####

#read data
mixture_all <- create_mixture_matrix(counts_filenames, mapping, data_file_list)
mixture_healthy <- create_mixture_matrix(counts_filenames_healthy, mapping, data_file_list_healthy)
mixture_unhealthy <- create_mixture_matrix(counts_filenames_unhealthy, mapping, data_file_list_unhealthy)


sample_number_list_healthy <- c("SRR12363244", "SRR12363245", "SRR12363246", 
                                "SRR12363247", "SRR12363248", "SRR16760002",
                                "SRR16760003", "SRR16760004", "SRR16760005",
                                "SRR16760006", "SRR16760007", "SRR16760008",
                                "SRR16760009", "SRR16760010", "SRR16760012")

sample_number_list_unhealthy <- c("SRR13632931", "SRR13632932", "SRR13632933",
                                  "SRR13632934", "SRR13632935", "SRR14788891",
                                  "SRR14788893", "SRR14788895", "SRR14788897",
                                  "SRR14788899" )

mixture_all_sort <- subset(mixture_all, select=c(SRR12363244,  SRR12363245 ,  SRR12363246 , 
                                                  SRR12363247 ,  SRR12363248 ,  SRR16760002 ,
                                                  SRR16760003 ,  SRR16760004 ,  SRR16760005 ,
                                                  SRR16760006 ,  SRR16760007 ,  SRR16760008 ,
                                                  SRR16760009 ,  SRR16760010 ,  SRR16760012 ,
                                                  SRR13632931 ,  SRR13632932 ,  SRR13632933 ,
                                                  SRR13632934 ,  SRR13632935 ,  SRR14788891 ,
                                                  SRR14788893 ,  SRR14788895 ,  SRR14788897 ,
                                                  SRR14788899 ))

colnamex_mix_all_sort <- data.frame(colnames(mixture_all_sort))

for (i in colnames(mixture_all_sort)) {
  if (i %in% sample_number_list_unhealthy) {
    colnamex_mix_all_sort$neu[colnamex_mix_all_sort$colnames.mixture_all_sort. ==i] <- paste0(i, "_U")
  }
  else
  {
    colnamex_mix_all_sort$neu[colnamex_mix_all_sort$colnames.mixture_all_sort. ==i] <- paste0(i, "_H")
  }
}

bab <- as.vector(colnamex_mix_all_sort$neu)
colnames(mixture_all_sort) <- bab

#### output table latex:
latextable_uh <- data.frame(Probenname = colnames(mixture_unhealthy),
                          Repository = rep("NCBI/BioProject", 10),
                          BioProject = rep(c("PRJNA649979","PRJNA736971"), each = 5),
                          ReadLength = rep(c(101, 150), each = 5),
                          seqdesign = rep(c("single", "paired"), each = 5))

rownames(latextable_uh) <- latextable_uh$Probenname

for (i in latextable_uh$Probenname) {
  latextable_uh[i, "Genes"] <- length(which(mixture_all[,i]!=0))
}

latextable_h <- data.frame(Probenname = colnames(mixture_healthy),
                            Repository = rep("NCBI/BioProject", 15),
                            BioProject = rep(c("PRJNA649979","PRJNA764684"), c(5,10)),
                            ReadLength = rep(c(101, 51), c(5,10)),
                            seqdesign = rep(c("single", "paired"), c(5,10)))

rownames(latextable_h) <- latextable_h$Probenname

for (i in latextable_h$Probenname) {
  latextable_h[i, "Genes"] <- length(which(mixture_all[,i]!=0))
}






####sc_signature_matrix:#####


#data read:
data_file_name <- "/home/eva/Documents/Alexander/endegelaende/test/c_single_cell_rna_data_hongmei.txt"
raw_data_read <- txt_to_dgCMatrix(data_file_name) 

data_file_name <- "/home/eva/Documents/Alexander/endegelaende/test/c_tsne_df_von_xls.csv"
spaltenbeschriftung <- colnames_from_txt(data_file_name)
reihenbeschriftung <- rownames_from_txt(data_file_name)
c_tsne_df <- as.data.frame(data_from_txt(data_file_name, spaltenbeschriftung, reihenbeschriftung))
c_tsne_df <- c_tsne_df[, -1]
c_tsne_df <- c_tsne_df[-1,]


raw_data_read_orig <- raw_data_read

for (j in 1:ncol(raw_data_read)) {
  colnames(raw_data_read)[j] <- c_tsne_df[colnames(raw_data_read)[j], "Apcluster"]
}

raw_data_read <- raw_data_read[,(colnames(raw_data_read) %in% c_tsne_df$Apcluster)]
#raw_data_read[1:4,1:4]
raw_data_missing <- raw_data_read[,!(colnames(raw_data_read) %in% c_tsne_df$Apcluster)]
#as.factor(raw_data_read)

#check geneexpression !=0

#raw_data_read <- raw_data_read[,colnames(raw_data_read) [colSums(raw_data_read)!=0]]
any(!(colSums(raw_data_read) > 0))
sc_for_signature_matrix <- raw_data_read[,!(colnames(raw_data_read) %in% c("EVT_24W_1", "EVT_24W_2"))]

raw_data_read[1:10,1:18]
raw_data_read_orig[1:10,1:18]
sc_for_signature_matrix[1:10,1:18]

####sc_signature_matrix_merged:#####
####beschriftungen:

merged_rules <- unique(data.frame(genes = c_tsne_df$Apcluster, class = c_tsne_df$Merge_Type))


merged_classes <- as.data.frame(colnames(raw_data_read))  
colnames(merged_classes) <- "genes"

merged_classes <- left_join(merged_classes, merged_rules)
sum(summary(as.factor(merged_classes$class)))
merged_classes$genes <- NULL

save_gene_subset_file_for_cs(merged_classes)

sc_for_signature_matrix_merged <- raw_data_read
colnames(sc_for_signature_matrix_merged) <- merged_classes$class

####nonsense matrix:####
werte_fuer_nonsense <- as.vector(raw_data_read)
set.seed(1)
sc_for_signature_matrix_merged_max_nonsense <- as.sparse(matrix(sample(werte_fuer_nonsense, length(werte_fuer_nonsense), replace = FALSE),
                                                                nrow = 20866, ncol = 1449, dimnames = list(
                                                                  rownames(sc_for_signature_matrix_merged),
                                                                  colnames(sc_for_signature_matrix_merged))))
nonsensecolnames <- rep_len(LETTERS[1:4], length(colnames(sc_for_signature_matrix_merged_max_nonsense)))
colnames(sc_for_signature_matrix_merged_max_nonsense) <- nonsensecolnames

rm(merged_rules, merged_classes, raw_data_missing, nonsensecolnames, reihenbeschriftung, 
   spaltenbeschriftung, werte_fuer_nonsense, data_file_list_unhealthy, data_file_list_healthy, data_file_list,
   data_file_name)

#####seurat####

###umap sc data inkl w24:###


fuer_pbmc <- raw_data_read
cbmc <- CreateSeuratObject(counts = fuer_pbmc,
                           # min.cells = 3, min.features = 200,
                           names.field = 1,
                           # names.delim = ".",
                           # meta.data = c_tsne_df,
                           row.names = rownames(fuer_pbmc))
set.seed(42)
# perform visualization and clustering steps
cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = TRUE)
#ElbowPlot(cbmc)
cbmc <- FindNeighbors(cbmc, dims = 1:50)
cbmc <- FindClusters(cbmc, resolution = 0.54, verbose = TRUE) #resolution: 0.025 fuer 7 cluster, 1.09 fuer 14 cluster, 0.724 fuer 1 bei nonsense
cbmc <- RunUMAP(cbmc, dims = 1:50)
cbmc <- RunTSNE(cbmc, dims = 1:50)

#umap:
cbmc_plot <- DimPlot(cbmc, label = TRUE)+
  ggtitle("UMAP plot - sig.matrix mit w24")
cbmc_plot

png(filename = paste0("umap_sc_","sig_mat_mit_W24",".png"))
cbmc_plot
dev.off()

#tsne:
cbmc[["umap"]] <- NULL
tsne_plot <- DimPlot(cbmc, label = TRUE)+
  ggtitle("t-SNE plot - signature matrix inkl. w24")
png(filename = paste0("tsne_sc_","sig_mat_mit_W24",".png"))
tsne_plot
dev.off()

# #PCA:
# cbmc[["tsne"]] <- NULL
# DimPlot(cbmc, label = TRUE)
# 

###umap fuer nonsense cluster:1###

fuer_pbmc <- sc_for_signature_matrix_merged_max_nonsense
cbmc <- CreateSeuratObject(counts = fuer_pbmc,
                           # min.cells = 3, min.features = 200,
                           names.field = 1,
                           # names.delim = ".",
                           # meta.data = c_tsne_df,
                           row.names = rownames(fuer_pbmc))
set.seed(42)
# perform visualization and clustering steps
cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = TRUE)
#ElbowPlot(cbmc)
cbmc <- FindNeighbors(cbmc, dims = 1:50)
cbmc <- FindClusters(cbmc, resolution = 0.54, verbose = TRUE) #resolution: 0.025 fuer 7 cluster, 1.09 fuer 14 cluster, 0.54 fuer 1 bei nonsense
cbmc <- RunUMAP(cbmc, dims = 1:50)
cbmc <- RunTSNE(cbmc, dims = 1:50)

#umap:
cbmc_plot <- DimPlot(cbmc, label = TRUE)+
  ggtitle("UMAP plot - nonsense")
cbmc_plot

png(filename = paste0("umap_sc_","nonsense",".png"))
cbmc_plot
dev.off()

#tsne:
cbmc[["umap"]] <- NULL
tsne_plot <- DimPlot(cbmc, label = TRUE)+
  ggtitle("t-SNE plot - nonsense")
png(filename = paste0("tsne_sc_","nonsense",".png"))
tsne_plot
dev.off()


# #PCA:
# cbmc[["tsne"]] <- NULL
# DimPlot(cbmc, label = TRUE)
# 

###umap fuer sig, mit einstellungen nonsense###
fuer_pbmc <- sc_for_signature_matrix
cbmc <- CreateSeuratObject(counts = fuer_pbmc,
                           # min.cells = 3, min.features = 200,
                           names.field = 1,
                           # names.delim = ".",
                           # meta.data = c_tsne_df,
                           row.names = rownames(fuer_pbmc))
set.seed(42)
# perform visualization and clustering steps
cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = TRUE)
#ElbowPlot(cbmc)
cbmc <- FindNeighbors(cbmc, dims = 1:50)
cbmc <- FindClusters(cbmc, resolution = 0.54, verbose = TRUE) #resolution: 0.025 fuer 7 cluster, 1.09 fuer 14 cluster, 0.724 fuer 1 bei nonsense
cbmc <- RunUMAP(cbmc, dims = 1:50)
cbmc <- RunTSNE(cbmc, dims = 1:50)



#umap:
cbmc_plot <- DimPlot(cbmc, label = TRUE)+
  ggtitle("UMAP plot - sig.matrix wie nonsense")
cbmc_plot

png(filename = paste0("umap_sc_","sig_mat_wie_nonsense",".png"))
cbmc_plot
dev.off()

#tsne:
cbmc[["umap"]] <- NULL
tsne_plot <- DimPlot(cbmc, label = TRUE)+
  ggtitle("t-SNE plot - signature matrix")
png(filename = paste0("tsne_sc_","sig_mat_wie_nonsense",".png"))
tsne_plot
dev.off()

# #PCA:
# cbmc[["tsne"]] <- NULL
# DimPlot(cbmc, label = TRUE)
# 


####daten zum pruefen der sig.matrix:####


#mixture_test_sample <- raw_data_read_orig[,sample(ncol(raw_data_read_orig), size = 25)]
mixture_test_sample <- sc_for_signature_matrix[,sample(ncol(sc_for_signature_matrix), size = 25)]
mixture_test_sample <- as.data.frame(mixture_test_sample)
mixture_test_sample_colnames <- colnames(mixture_test_sample)
for (i in 1:length(mixture_test_sample_colnames)) {
  mixture_test_sample_colnames[i] <- gsub("_","", paste0(mixture_test_sample_colnames[i],"sample",i))
}
colnames(mixture_test_sample) <- mixture_test_sample_colnames
mixture_test_sample$tx.gene_name <- rownames(mixture_test_sample)

testbasis <- mapping
testbasis$tx.tx_id <- NULL
testbasis <- unique(testbasis)
rownames(testbasis) <- testbasis$tx.gene_name

mixture_test_like_all <- full_join(testbasis, mixture_test_sample)
mixture_test_like_all[is.na(mixture_test_like_all)] <- 0
rownames(mixture_test_like_all) <- mixture_test_like_all$tx.gene_name
mixture_test_like_all$tx.gene_name <- NULL
mixture_test_sample_aus_sig_matrix <- subset(mixture_test_sample, select = -tx.gene_name)


summary(as.factor(rownames(mixture_test_like_all) %in% rownames(sc_for_signature_matrix)))




#roh-daten in .txt für cibersortx speichern:
save_mixture_file_for_cs(mixture_all)
save_mixture_file_for_cs(mixture_healthy)
save_mixture_file_for_cs(mixture_unhealthy)
save_mixture_file_for_cs(mixture_test_sample_aus_sig_matrix)
save_mixture_file_for_cs(mixture_test_like_all)
save_mixture_file_for_cs(mixture_all_sort)


save_mixture_file_for_cs(sc_for_signature_matrix_merged_max_nonsense)
save_mixture_file_for_cs(sc_for_signature_matrix_merged)
save_mixture_file_for_cs(sc_for_signature_matrix)
# 
# write.table(latextable_h, file = "latextable_h.txt", sep = " & ", quote = FALSE,
#             row.names = FALSE, eol = "\\ \r")
# write.table(latextable_uh, file = "latextable_uh.txt", sep = " & ", quote = FALSE,
#             row.names = FALSE, eol = "\\ \r")



rm(list = ls())
gc()


