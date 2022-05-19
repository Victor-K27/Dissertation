
#loading libraries
library(tibble)
library(qqman)

#setting directory and reading the files
setwd("C:/Users/victo/Desktop/Fourth year project/Data/141 genes eQTL analysis")
downregulated_genes <- read.delim("eqtls_downregulated_markers.corrected.txt")
upregulated_genes <- read.delim("eqtls_upregulated_markers.corrected.txt")

#adding downregulated and upregulated status for each gene in the dataframe
downregulated_genes <- add_column(downregulated_genes, "downregulated", .after = "gene_name")
colnames(downregulated_genes)[4] <- "gene_status"
upregulated_genes <- add_column(upregulated_genes, "upregulated", .after = "gene_name")
colnames(upregulated_genes)[4] <- "gene_status"

#extracting vectors with names of the genes for heatmap on GTEX portal
downregulated_genes_names <- unique(downregulated_genes$gene_name)
upregulated_genes_names <- unique(upregulated_genes$gene_name)
write(downregulated_genes_names, "downregulated_genes_names.txt")
write(upregulated_genes_names, "upregulated_genes_names.txt")

    #Genetic control analysis- "Genetic control.R" 
#combine the dataframes
all_genes_dataframe <- rbind(downregulated_genes,upregulated_genes)

#making a list of all the genes
myfiles <- split.data.frame(all_genes_dataframe, all_genes_dataframe$gene_name)

#making a vector with all unique eQTL (significant or not)
all_eQTL <- unique(all_genes_dataframe$rs_id_dbSNP151_GRCh38p7)

#making a vector with unique significant eQTL
uniqueeQTL <- unique(all_genes_dataframe[all_genes_dataframe$qval < 0.05,]$rs_id_dbSNP151_GRCh38p7)

#compare eQTLs between files
eQTLdataframe <- data.frame(uniqueeQTL)
for (i in 1:length(myfiles)){
  eQTLdataframe[myfiles[[i]]$gene_name[1]] <- uniqueeQTL %in% myfiles[[i]]$rs_id_dbSNP151_GRCh38p7
}

#naming the rows with the name of the unique eQTLs
rownames(eQTLdataframe)= eQTLdataframe$uniqueeQTL
#removing the first column so that we have gene 1 in column 1
eQTLdataframe= eQTLdataframe[-c(1)]

#adding a column with the sum of the genes in which the eQTL is found at (e.g. if eQTL x is found in genes a,b and c, the sum would be 3)
eQTLdataframe["sum"] <- rep(0,nrow(eQTLdataframe))
for (i in 1:nrow(eQTLdataframe)){
  eQTLdataframe[i, "sum"] <- sum(eQTLdataframe[i,])
}

#adding a column with the eQTL IDs in the dataframe
eQTLdataframe["eQTL IDs"] <- rownames(eQTLdataframe)
#isolating the names of the eQTLs that control more than one gene
eQTL_names <- eQTLdataframe[eQTLdataframe$sum > 1,]$`eQTL IDs`
#making a dataframe with data for the significant eQTLs
significant_eQTLdataframe <- all_genes_dataframe[all_genes_dataframe$rs_id_dbSNP151_GRCh38p7 %in% eQTL_names,]
#remove the "." 
significant_eQTLdataframe <- significant_eQTLdataframe[!significant_eQTLdataframe$rs_id_dbSNP151_GRCh38p7 == ".",]
#select qval < 0.05
significant_eQTLdataframe <- significant_eQTLdataframe[which(significant_eQTLdataframe$qval < 0.05),]
eQTL_two_genes_each <- unique(significant_eQTLdataframe$rs_id_dbSNP151_GRCh38p7)
write(eQTL_two_genes_each, "eQTLs_two_genes_each.txt")

#writing a file with the significant_eQTLs in the 141 genes
if (askYesNo("Do you want to save significant_eQTLdataframe?", )==TRUE){
  write.csv2(significant_eQTLdataframe,"/Users/victorkonstantellos/Desktop/Fourth year project/Data/135 genes analysis/significant_eQTLdataframe.csv", row.names = TRUE)
}

      #Manhattan plot for 141 genes
#making the chr_integer column out of chr 
all_genes_dataframe$chr_integer <- as.integer(gsub("[^[:digit:]]", "", all_genes_dataframe$chr))

#creating a dataframe for manhattan plot
manhattan_dataframe <- data.frame(SNP = all_genes_dataframe$rs_id_dbSNP151_GRCh38p7, CHR = all_genes_dataframe$chr_integer, 
                                  BP = all_genes_dataframe$variant_pos, P = all_genes_dataframe$pval_true_df)

#keeping only the unique rows in the dataframe
manhattan_dataframe <- na.omit(manhattan_dataframe[!duplicated(manhattan_dataframe$SNP),])

#plotting the Manhattan plot
manhattan(x= manhattan_dataframe, main = "Manhattan Plot, Unique eQTLs for 141 genes across 49 Tissues", 
          genomewideline = -log(0.05), suggestiveline = FALSE, highlight = eQTL_two_genes_each)

