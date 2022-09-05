# library(readxl)
 library(Seurat)
 library(dplyr)
 library(ggplot2)
# library(tidyr)
# library(ggridges)
# library(readxl)
# library(DESeq2)
# library(SeuratObject)
# library(patchwork)
# library(umap)
# library(reticulate)
# library(tidyr)
# library(stats)
# library(forcats)
# ##library(hrbrthemes),
# library(viridis)
 library(gplots)



####functions####

# to plot graphics into file:
plot_to_pdf <- function(plot_titel){
  pdf(file = plot_titel, width = 34, height = 12)
  dev.off()
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

significantly_differentially_expressed_genes_between_groups <- function(geps1, geps2, stderr1, stderr2){
  for (i in 1:ncol(geps1)){
    vBetaZ <- sapply(1:nrow(geps1), function(j) (geps1[j,i]-geps2[j,i])/sqrt(stderr1[j,i]^2+stderr2[j,i]^2))
    ZPs <- 2*pnorm(-abs(vBetaZ))
    Zqvals <- p.adjust(ZPs, method = "BH")
  }
  df <- data.frame(gene_name= rownames(geps1), Pvalue = Zqvals)
  return(Zqvals)
}



####directory,files, etc####



#setting workingdirectory 
### for manual run:
setwd("~/Documents/Alexander/endegelaende/test/cs_in")
#setwd("~/Downloads/")

#setwd(get_script_path())

filename_data_geps_26 <- "CIBERSORTx_Job26_output/CIBERSORTxGEP_Job26_GEPs.txt"
filename_data_geps_27<- "CIBERSORTx_Job27_output/CIBERSORTxGEP_Job27_GEPs.txt"
filename_data_geps_26_filtered <- "CIBERSORTx_Job26_output/CIBERSORTxGEP_Job26_GEPs_Filtered.txt"
filename_data_geps_27_filtered <- "CIBERSORTx_Job27_output/CIBERSORTxGEP_Job27_GEPs_Filtered.txt"
filename_data_stderr_26 <- "CIBERSORTx_Job26_output/CIBERSORTxGEP_Job26_GEPs_StdErrs.txt"
filename_data_stderr_27 <- "CIBERSORTx_Job27_output/CIBERSORTxGEP_Job27_GEPs_StdErrs.txt"
#filename_data_geps_18 <- "CIBERSORTx_Job18_output/CIBERSORTxGEP_Job18_GEPs.txt"
#filename_data_geps_18_filtered <- "CIBERSORTx_Job18_output/CIBERSORTxGEP_Job18_GEPs_Filtered.txt"

geps_1 <- as.matrix(txt_to_dgCMatrix(filename_data_geps_26))
geps_2 <- as.matrix(txt_to_dgCMatrix(filename_data_geps_27))
geps_1_filtered <- as.matrix(txt_to_dgCMatrix(filename_data_geps_26_filtered))
geps_2_filtered <- as.matrix(txt_to_dgCMatrix(filename_data_geps_27_filtered))
stderr_1 <- as.matrix(txt_to_dgCMatrix(filename_data_stderr_26))
stderr_2 <- as.matrix(txt_to_dgCMatrix(filename_data_stderr_27))
#geps_18 <- as.matrix(txt_to_dgCMatrix(filename_data_geps_18))
#geps_18_filtered <- as.matrix(txt_to_dgCMatrix(filename_data_geps_18_filtered))



unfiltered <- significantly_differentially_expressed_genes_between_groups(geps_1, geps_2, stderr_1, stderr_2)
filtered <- significantly_differentially_expressed_genes_between_groups(geps_1_filtered, geps_2_filtered, stderr_1, stderr_2)  


###find interesting genes for hi-res ####
#from list:
# data_file_name <- "gene_list.txt"
# spaltenbeschriftung <- colnames_from_txt(data_file_name)
# reihenbeschriftung <- rownames_from_txt(data_file_name)
# gene_subset_fi <- as.data.frame(data_from_txt(data_file_name, spaltenbeschriftung, reihenbeschriftung))
# 
# #die gewünschten Gene einsetzen!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# set.seed(1)
# gene_subset_file <- sample(gene_subset_fi$ACE, 150, replace=FALSE)
# 
# save_gene_subset_file_for_cs(gene_subset_file)
# 


#from GEP most expressed:
# find_interesting_genes <- geps_18_filtered
# find_interesting_genes <- as.data.frame(find_interesting_genes)
# find_interesting_genes_stb <- subset (find_interesting_genes, find_interesting_genes$STB > 5, select = c(STB))#, CTB_8W_3,CTB_8W_2,CTB_8W_1))
# find_interesting_genes_stb$genename <-  rownames(find_interesting_genes_stb)
# find_interesting_genes <- geps_18_filtered
# find_interesting_genes <- as.data.frame(find_interesting_genes)
# find_interesting_genes_ctb <- subset (find_interesting_genes, find_interesting_genes$CTB > 170, select = c(CTB))#, CTB_8W_3,CTB_8W_2,CTB_8W_1))
# find_interesting_genes_ctb$genename <-  rownames(find_interesting_genes_ctb)
# 
# #find_interesting_genes[is.na(find_interesting_genes)] <- 0
# #find_interesting_genes <- as.data.frame(cbind(find_interesting_genes, summe_exp = rowSums(find_interesting_genes)))
# #find_interesting_genes_2 <- find_interesting_genes[with(find_interesting_genes, order(-summe_exp)),]
# #find_interesting_genes_2 <- find_interesting_genes_2[1:1000,]
# #set.seed(1)
# #gene_subset_file_stb_ctb <- sample(rownames(find_interesting_genes_2), 1000, replace=FALSE)
# gene_subset_file_stb_ctb <- full_join(find_interesting_genes_ctb, find_interesting_genes_stb, )
# gene_subset_file_stb_ctb <- subset(gene_subset_file_stb_ctb, select = c(genename))
# 
# data_file_name <- "/home/eva/Downloads/gene_list.txt"
# spaltenbeschriftung <- colnames_from_txt(data_file_name)
# reihenbeschriftung <- rownames_from_txt(data_file_name)
# gene_subset_fi <- as.data.frame(data_from_txt(data_file_name, spaltenbeschriftung, reihenbeschriftung))
# colnames(gene_subset_fi) <- c("genename")
# gene_subset_file_stb <- full_join(gene_subset_fi, find_interesting_genes_stb )
# gene_subset_file_stb <- gene_subset_file_stb$genename
# 
# save_gene_subset_file_for_cs(gene_subset_file_stb)


## find significant genes between gep1 & 2####
df <- data.frame(gene_name= rownames(geps_1), Pvalue_unfiltered = unfiltered, Pvalue_filtered = filtered)
significantly_deg <- df
rownames(significantly_deg) <- significantly_deg$gene_name
significantly_deg[is.na(significantly_deg)] <- 1
#any(is.na(significantly_deg))
significantly_deg_all <- significantly_deg %>%
  mutate_if(is.numeric,
            round,
            digits = 5)


significantly_deg <- subset(significantly_deg_all, significantly_deg_all$Pvalue_filtered < 0.5 | significantly_deg_all$Pvalue_unfiltered < 0.5)

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

heatmapSigGenes <- (plot_table + scale_fill_viridis_c(option = "D") + labs(x = "Filtered/unfiltered", y = NULL))

png(file = "GroupGepHeatmap.png", width = 640, height = 640)
heatmapSigGenes
dev.off()


#gene zwischen gesund und krank auf einer tSNE/PCA darstellen?
#pdf(file = "GEP_heatmap_from_r",width = 24, height = 24)



genelist_for_cs <- subset(significantly_deg_all, significantly_deg_all$Pvalue_filtered < 1 | significantly_deg_all$Pvalue_unfiltered < .92)

save_gene_subset_file_for_cs(genelist_for_cs)

save_mixture_file_for_cs(significantly_deg)

# 
# pdf(file = "fraction_plots.pdf",width = 14, height = 14)
# boxplot_celltypes
# density_binline_celltypes
# density_celltypes
# dev.off()

rm(list = ls())
gc()
