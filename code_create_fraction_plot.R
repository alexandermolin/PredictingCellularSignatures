
library(readxl)
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggridges)
library(pheatmap)
library(tidyverse)
library(ggprism)
library(patchwork)
library(magrittr)


## ---- funktionen


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



## ---- directories



#setting workingdirectory 
### for manual run:
setwd("~/Documents/Alexander/endegelaende/test/cs_in/mixtureAllsignatureMatrix/")
#setwd("~/Documents/Alexander/endegelaende/test/cs_in/mixtureAllsignatureNonsense/")
#setwd("~/Documents/Alexander/endegelaende/test/cs_in/")
#setwd("~/Downloads/")


#setwd(get_script_path())


#list of files
data_file_list <- list.files(all.files = FALSE, full.names = FALSE, pattern = "^CIBERSORTx_Job")
data_file_list_names <- gsub("CIBERSORTx_","",data_file_list)
data_file_list_names <- gsub("_Adjusted.txt","",data_file_list_names)

fraction_daten_list <- list()

for (i in 1:length(data_file_list)) {
  
  spaltenbeschriftung <- colnames_from_txt(data_file_list[i])
  reihenbeschriftung <- rownames_from_txt(data_file_list[i])
  raw_data_read <- as.data.frame(data_from_txt(data_file_list[i], spaltenbeschriftung, reihenbeschriftung))
  raw_data_read <- raw_data_read[, -1]
  raw_data_read <- raw_data_read[-1,]
  
fraction_daten_list[[i]] <- raw_data_read
names(fraction_daten_list)[i] <- data_file_list_names[i]
}





# yea <-  7
# 
# ##change for nonsense:!!!!!!!!!!!!!!!!!!!!!1
# filename_data_fractions_adjusted <- paste0("CIBERSORTx_Job",yea ,"_Adjusted.txt")
# #filename_data_fractions_adjusted <- paste0("CIBERSORTx_Job",yea ,"_Results.txt")
# 
# 
# cs_adjusted <- as.data.frame(txt_to_dgCMatrix(filename_data_fractions_adjusted))


sample_number_list_healthy <- c("SRR12363244", "SRR12363245", "SRR12363246", 
                            "SRR12363247", "SRR12363248", "SRR16760002",
                            "SRR16760003", "SRR16760004", "SRR16760005",
                            "SRR16760006", "SRR16760007", "SRR16760008",
                            "SRR16760009", "SRR16760010", "SRR16760012")

sample_number_list_unhealthy <- c("SRR13632931", "SRR13632932", "SRR13632933",
                              "SRR13632934", "SRR13632935", "SRR14788891",
                              "SRR14788893", "SRR14788895", "SRR14788897",
                              "SRR14788899" )


zellen_reihenfolge_plot_merged <- c("EVT_8W",
                                    "STB", 
                                    "CTB",
                                    "Marc", 
                                    "Blood",
                                    "Mes",
                                    "P.value", "Correlation","RMSE")

zellen_reihenfolge_plot <- c("EVT_8W_3","EVT_8W_2","EVT_8W_1",
                             "STB_8W", 
                             "CTB_8W_3", "CTB_8W_2", "CTB_8W_1",
                             "Macro_2", "Macro_1", 
                             "blood.cell",
                             "Mes_2","Mes_1",
                             "P.value", "Correlation","RMSE")

zellen_reihenfolge_plot_nonsense <- c("A", "B", "C", "D",
                                      "P.value", "Correlation","RMSE")

###change if merged!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!:
  
zellen_reihenfolge_fuer_plot <- zellen_reihenfolge_plot



###plot####

j <- 1


celltypes_for_violin <- as.data.frame(fraction_daten_list[j])
colnames(celltypes_for_violin) <- gsub(paste0(data_file_list_names[j],"."),"", colnames(celltypes_for_violin))
samplenames <- rownames(celltypes_for_violin)
celltypes_for_violin <- as.data.frame(sapply(celltypes_for_violin, as.numeric))
celltypes_for_violin$samplenames <- samplenames
rownames(celltypes_for_violin) <- samplenames
celltypes_for_violin$Health <- c(rep("healthy",length(celltypes_for_violin[,1])))

###nicht bie testing:
for (i in rownames(celltypes_for_violin)) {
  if (i %in% sample_number_list_unhealthy) {
    celltypes_for_violin$Health[rownames(celltypes_for_violin)==i] <- "unhealthy"
  }
  }

for (i in rownames(celltypes_for_violin)) {
  if (i %in% sample_number_list_unhealthy) {
    celltypes_for_violin$samplenames[rownames(celltypes_for_violin)==i] <- paste0(i, "_U")
  }
    else
    {
      celltypes_for_violin$samplenames[rownames(celltypes_for_violin)==i] <- paste0(i, "_H")
  }
}

#celltypes_for_violin$samplenames <- rownames(celltypes_for_violin)
#celltypes_for_violin$Pvalue <- celltypes_for_violin$`P-value`
#celltypes_for_violin$`P-value` <- NULL
data_pivot <- as_tibble(subset(celltypes_for_violin, select = c(-P.value, -Correlation, -RMSE)))
data_pivot_2 <- pivot_longer(data_pivot, -c(samplenames, Health), names_to = "Celltypes", values_to = "Fraction")
#data_pivot_3 <- as_tibble(subset(celltypes_for_violin, select = c(-Health, -samplenames)))
data_pivot_2$Celltypes <- factor(data_pivot_2$Celltypes, ordered = T, levels = zellen_reihenfolge_fuer_plot)

data_pivot_3 <- pivot_longer(celltypes_for_violin, -c(Health, samplenames),names_to = "Celltypes", values_to = "Fraction")

# Most basic violin chart

boxplotCelltypes <- ggplot(data_pivot_2, aes(x = Celltypes, y = Fraction, fill = Health)) +
  geom_boxplot(outlier.shape = 1) +
  geom_point(shape = 4, size = 0.4, alpha = 0.7) # width = 0.25,

boxplotCelltypes

densityBinLine <- ggplot(data_pivot_2, aes(y = Celltypes, x = Fraction, fill = Health)) +
  #geom_boxplot(outlier.shape = 1) +
  #geom_point(shape = 4, width = 0.25, size = 0.4, alpha = 0.7)
  geom_density_ridges(alpha = 0.4, stat = "binline", bins = 20)


densityCelltype <- ggplot(data_pivot_2, aes(y = Celltypes, x = Fraction, fill = Health)) +
  #geom_boxplot(outlier.shape = 1) +
  #geom_point(shape = 4, width = 0.25, size = 0.4, alpha = 0.7)
  geom_density_ridges(alpha = 0.4)



# data_for_plot_table <- as.data.frame(celltypes_for_violin)
# data_for_plot_table$Health <- as.factor(data_for_plot_table$Health)
# data_for_plot_table <- data_for_plot_table[,]
data_pivot_3 <- data_pivot_3 %>%
  mutate_if(is.numeric,
            round,
            digits = 3)



data_pivot_3$Celltypes <- factor(data_pivot_3$Celltypes, ordered = T, levels = zellen_reihenfolge_fuer_plot)

###############nicht bei matrix testing:
sample_number_list_healthy_new <- paste(sample_number_list_healthy, "_H", sep = "") 
sample_number_list_unhealthy_new <- paste(sample_number_list_unhealthy, "_U", sep = "") 

data_pivot_3$samplenames <- factor(data_pivot_3$samplenames, ordered = T,levels = c(sample_number_list_healthy_new,
                                                                                    sample_number_list_unhealthy_new))

plot_table <- ggplot(data_pivot_3, aes(x=Celltypes, y=samplenames, fill=Fraction))+
  geom_tile(color = "black", lwd = 1.5, linetype = 1) + 
  geom_text(aes(label = Fraction), color = "white", size = 8) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20, title = "Fractions")) +
  ylim(rev(levels(as.factor(data_pivot_3$samplenames)))) + 
  xlim(levels(as.factor(data_pivot_3$Celltypes)))

fractionsTable <- (plot_table + scale_fill_viridis_c(option = "D") +
                     labs(x = "Celltypes", y = NULL) +
                     theme(axis.text = element_text(size = 20, angle=0),
                           axis.title = element_text(size = 22, angle=0)))

png(file = paste0("fraction", data_file_list_names[j], "Table.png"),width = 2000, height = 960)
fractionsTable
dev.off()


####for testing the matrix:

celltypes_for_violin[,"max"] <- apply(celltypes_for_violin[1:17], 1, max)

celltypes_for_violin$max <- do.call(pmax, celltypes_for_violin[1:12])

celltypes_for_violin$max <- celltypes_for_violin[1:12][cbind(seq_len(nrow(celltypes_for_violin)), max.col(celltypes_for_violin[1:12]))]
celltypes_for_violin$max_cell <- colnames(celltypes_for_violin)[max.col(celltypes_for_violin[1:12])]
celltypes_for_violin$richtig <- rep(1,25)
celltypes_for_violin[c(17,18,19,25),"richtig"] <- 0
data_pivot_5 <- pivot_longer(celltypes_for_violin, -c(Health, samplenames, max_cell),names_to = "Celltypes", values_to = "Fraction")
data_pivot_5 <- data_pivot_5 %>%
  mutate_if(is.numeric,
            round,
            digits = 3)

zellen_reihenfolge_fuer_plot_testing <- c("EVT_8W_3","EVT_8W_2","EVT_8W_1",
                             "STB_8W", 
                             "CTB_8W_3", "CTB_8W_2", "CTB_8W_1",
                             "Macro_2", "Macro_1", 
                             "blood.cell",
                             "Mes_2","Mes_1",
                             "P.value", "Correlation","RMSE", "max", "richtig")

data_pivot_5$Celltypes <- factor(data_pivot_5$Celltypes, ordered = T, levels = zellen_reihenfolge_fuer_plot_testing)
data_pivot_5$samplenames <- factor(data_pivot_5$samplenames)


plot_table <- ggplot(data_pivot_5, aes(x=Celltypes, y=samplenames, fill=Fraction))+
  geom_tile(color = "black", lwd = 1.5, linetype = 1) + 
  geom_text(aes(label = Fraction), color = "white", size = 8) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20, title = "Anteil")) +
  ylim(rev(levels(as.factor(data_pivot_5$samplenames)))) + xlim(levels(as.factor(data_pivot_5$Celltypes)))

fractionsTable <- (plot_table + scale_fill_viridis_c(option = "D") + 
                     labs(x = "Celltypes", y = NULL) +
                     theme(axis.text = element_text(size = 18, angle=0),
                           axis.title = element_text(size = 22, angle=0)))

png(file = paste0("fraction", data_file_list_names[j], "Table.png"),width = 2000, height = 960)
fractionsTable
dev.off()


#### ende testing


###two sample wilcoxon-test:####

data_wilcoxon <- celltypes_for_violin[,c(1:12,17)]
data_wilcoxon$Health <- as.factor(data_wilcoxon$Health)
rownames(data_wilcoxon) <- NULL
results_wilcoxon <- list()

for (w in 1:12) {
  x <- data_wilcoxon[[w]]
  y <- data_wilcoxon$Health
  results_wilcoxon[[w]] <- wilcox.test(x~y, exact = TRUE, conf.int = TRUE)
}

a <- unlist(cbind(results_wilcoxon[[1]][3],
       results_wilcoxon[[2]][3],
       results_wilcoxon[[3]][3],
       results_wilcoxon[[4]][3],
       results_wilcoxon[[5]][3],
       results_wilcoxon[[6]][3],
       results_wilcoxon[[7]][3],
       results_wilcoxon[[8]][3],
       results_wilcoxon[[9]][3],
       results_wilcoxon[[10]][3],
       results_wilcoxon[[11]][3],
       results_wilcoxon[[12]][3]))

names(a) <- colnames(celltypes_for_violin[1:12])
a <- round(a, digits = 5)
a

label_dataf <- data.frame(celltype = names(a),
                          pval=a)

dfPval <- data.frame(cellt=zellen_reihenfolge_plot[1:12])
pvalfurplot <- c(0.13159,0.47773,0.85968,0.00090,0.84896,0.71791,
                 0.02899,0.97783,0.00141,0.75776,0.00004,0.00692)




boxplotCelltypes <- ggplot(data_pivot_2, aes(x = Celltypes, y = Fraction, fill = Health)) +
  geom_boxplot(outlier.shape = 1) +
  geom_point(shape = 4, size = 0.4, alpha = 0.7) +
  annotate(geom = "text", x= 0.8, y= 0.95,size=8, label= "p-values:")+
  annotate(geom = "text", x= 1:12, y= 0.90,size=8, label=c("0.13159","0.47773","0.85968","0.00090",
                                                   " 0.84896","0.71791","0.02899","0.97783",
                                                   "0.00141","0.75776","0.00004","0.00692"))+
  theme(axis.text = element_text(size = 18, angle=0),
        axis.title = element_text(size = 22, angle=0))



png(file = paste0("fraction", data_file_list_names[j], "Boxplotfaa.png"),width = 2000, height = 960)
boxplotCelltypes
dev.off()

boxplotCelltypes 

# str(results_wilcoxon)
# data_wilcoxon[[1]]
# 
# names(results_wilcoxon) <- colnames(celltypes_for_violin[1:12])
# 
# result_wilcoxon[1]
# 
# out_res_wilcoxon <- lapply(1:12, function(x) 
#   pairwise.wilcox.test(data_wilcoxon[[x]], data_wilcoxon$Health, p.adjust.method = "bonf"))
# names(out_res_wilcoxon) <- names(data_wilcoxon)[1:12]
# out_res_wilcoxon
# sapply(out_res_wilcoxon, function(x) {
#   p <- x$p.value
#   n <- outer(rownames(p), colnames(p), paste, sep='v')
#   p <- as.vector(p)
#   names(p) <- n
#   p
# })

####P-Value vergleich table:#####


jobs <- c("Job7", "Job8", "Job9", "Job10", "Job11")
permutations <- c(0, 50, 100, 500, 1000)

rm(dfpvalue, dfpvalue2, data_pivot_4)

dfpvalue <- data.frame(fraction_daten_list[[1]]$`P-value`,
                          fraction_daten_list[[2]]$`P-value`,
                          fraction_daten_list[[3]]$`P-value`,
                          fraction_daten_list[[4]]$`P-value`,
                          fraction_daten_list[[5]]$`P-value`)
dfpvalue <- data.frame(sapply(dfpvalue, as.numeric))
colnames(dfpvalue) <- data_file_list_names
rownames(dfpvalue) <- rownames(fraction_daten_list[[1]])


dfpvalue_samplenames <- rownames(dfpvalue)

#data_pivot_3 <- pivot_longer(dfpvalue, cols = c(Job7,Job8, Job9, Job10, Job11),  names_to = "Job", values_to = "values")
#data_pivot_3$values <- as.numeric(data_pivot_3$values)
dfpvalue <- dfpvalue %>%
  mutate_if(is.numeric,
            round,
            digits = 5)
dfpvalue <- data.frame(Job7 = dfpvalue$Job7,
                       Job8 = dfpvalue$Job8,
                       Job9 = dfpvalue$Job9,
                       Job10 = dfpvalue$Job10,
                       Job11 = dfpvalue$Job11)
rownames(dfpvalue) <- dfpvalue_samplenames                       
colnames(dfpvalue) <- permutations

dfpvalue2 <- dfpvalue[rowSums(dfpvalue)>0,]
dfpvalue2$samplenames <- rownames(dfpvalue2)


data_pivot_4 <- pivot_longer(dfpvalue2, 
                             cols = as.character(permutations),
                             names_to = "Job", values_to = "values")
data_pivot_4$Job <- factor(data_pivot_4$Job, ordered = T,
                                 levels = permutations)


plot_table <- ggplot(data_pivot_4, aes(x=Job, y=samplenames, fill=values))+
  geom_tile(color = "black", lwd = 1.5, linetype = 1) +
  geom_text(aes(label = values), color = "white", size = 8) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20, title = "P-Values")) +
  ylim(rev(levels(as.factor(data_pivot_4$samplenames)))) + xlim(levels(as.factor(data_pivot_4$Job)))


fractionsTable <- (plot_table + scale_fill_viridis_c(option = "D") +
                     labs(x = "Permutations", y = "Samples mit P-Value < 0") +
                     theme(axis.text = element_text(size = 18, angle=0),
                           axis.title = element_text(size = 22, angle=0)))
                   
png(file = paste0("table_pValue_nonsense.png"),width = 1000, height = 460)
fractionsTable
dev.off()



rm(list = ls())
gc()

