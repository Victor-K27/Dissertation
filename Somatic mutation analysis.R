                          #run the 135 gene analysis before this code if you want to test which of the significant eQTLs match with somatic SNPs from the TCGA
#loading the packages                       
library(TCGAbiolinks)
library(maftools)
library(dplyr)
library(ggplot2)
library(readxl)
library(RTCGAToolbox)
library(rstatix)
library(ggpubr)

#loading the data
load("myEnvironment.RData")

#reading the Pancancer file
setwd("~/Desktop/Fourth year project/Data/Maftools")
pancancer_quiescence_groups <- read.csv2("PancancerQuiescenceGroupsCorrected.csv")

#loading the files from TCGA database- takes about 30 minutes
if (askYesNo("Do you want to load TCGA datasets") == TRUE){
  cancer_types <- c("LAML","ACC","BLCA","LGG","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","DLBC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THYM","THCA","UCS","UCEC","UVM")
  cancer_types_dataframe <- data.frame()
  for (i in cancer_types){
  maf <- GDCquery_Maf(i, pipelines = "mutect2") %>% read.maf
  df <- data.frame(maf@data)
  df$CancerStudy <- i
  cancer_types_dataframe <- rbind(cancer_types_dataframe, df)
  }
}
#"LCML" "CNTL" "FPPP" "MISC" did not have data

      ## eQTL analysis

#reading the cancer_type_dataframe: this has the short barcode as well. 
if (askYesNo("Do you want to read the cancer_types_dataframe from the TCGA_data?") == TRUE){
  cancer_types_dataframe <- read.csv2("TCGA_data.csv")
}

#creating a dataframe with the TCGA data that matched the SNP dataframe
TCGA_SNPs_merged <- merge(cancer_types_dataframe, significant_SNPdataframe, by.x = c("Chromosome", "Start_Position"), by.y = c("chr", "variant_pos"))
any(significant_SNPdataframe$variant_pos %in% cancer_types_dataframe$Start_Position) == TRUE

#creating a dataframe with the TCGA data that matched the all_genes_dataframe
TCGA_allgenes_merged <- merge(cancer_types_dataframe, all_genes_dataframe, by.x = c("Chromosome", "Start_Position"), by.y = c("chr", "variant_pos"))
any(all_genes_dataframe$variant_pos %in% cancer_types_dataframe$Start_Position) == TRUE


    ## Mutations analysis
#merging the cancer types dataframe with pancancer quiescne groups
cancer_types_dataframe_merged_pancancer <- merge(cancer_types_dataframe, pancancer_quiescence_groups, by = "Barcode2")

#splitting the samples into quiescent and fast cycling
quiescent_pancancer <- cancer_types_dataframe_merged_pancancer[cancer_types_dataframe_merged_pancancer$Group == "Highly Quiescent",]
fast.cycling_pancancer <- cancer_types_dataframe_merged_pancancer[cancer_types_dataframe_merged_pancancer$Group == "Fast Cycling",]
#reading them as maf files
quiescent_pancancer_maf <- read.maf(quiescent_pancancer)
fast.cycling_pancancer_maf <- read.maf(fast.cycling_pancancer)

coBarplot(m1 = quiescent_pancancer_maf, m2= fast.cycling_pancancer_maf, m1Name = "Highly Quiescent, Pancancer", m2Name = "Fast Cycling, Pancancer")

fisher_test_pancancer <- mafCompare(m1 = quiescent_pancancer_maf, m2 = fast.cycling_pancancer_maf, 
                                    m1Name = "Highly Quiescent", m2Name = "Fast Cycling", minMut = 5)
fisher_test_pancancer_results <- as.data.frame(fisher_test_pancancer)
fisher_test_pancancer_results$results.pval.adj <- p.adjust(fisher_test_pancancer_results$results.pval, method = "fdr")
forestPlot(fisher_test_pancancer, fdr = 5e-22, titleSize = 1.5, geneFontSize = 0.9)

#cancer-specific fisher test applied to all cancers
setwd("~/Desktop/Fourth year project/Data/Maftools/GDCdata_SNV")
cancer_files <- list.files(path = "~/Desktop/Fourth year project/Data/Maftools/GDCdata_SNV", pattern = ".gz")
sample_size_dataframe <- data.frame(Cancer= NA, N_quiescent= NA, N_fast.cycling= NA)
fisher_test_results <- list()
cancers_no_sig_genes <- c()
cancers_not_enough_samples <- c()
cancers_no_barcode_match <- c()
for (i in 1:3){
  data <- read.maf(cancer_files[i])
  data <- data@data
  
  data$Barcode2 <- sapply(as.character(data$Tumor_Sample_Barcode), function(x) paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
  
  if (any(data$Barcode2 %in% pancancer_quiescence_groups$Barcode2) == TRUE){
    data_merged_pancancer <- merge(data, pancancer_quiescence_groups, by = "Barcode2")
    
    data_quiescent_temporary <- data_merged_pancancer[data_merged_pancancer$Group == "Highly Quiescent",]
    data_fast.cycling_temporary <- data_merged_pancancer[data_merged_pancancer$Group == "Fast Cycling",]
    
    sample_size_dataframe <- rbind(sample_size_dataframe, c(Cancer= as.character(data_merged_pancancer$CancerType)[1], N_quiescent= length(unique(data_quiescent_temporary[data_quiescent_temporary$Group == "Highly Quiescent"]$Barcode2)), 
                                                            N_fast.cycling= length(unique(data_fast.cycling_temporary[data_fast.cycling_temporary$Group == "Fast Cycling",]$Barcode2))))
    
    if (nrow(data_quiescent_temporary) > 0 & nrow(data_fast.cycling_temporary) > 0){
      data_quiescent_temporary_maf <- read.maf(data_quiescent_temporary)
      data_fast.cycling_temporary_maf <- read.maf(data_fast.cycling_temporary)
      
      
    fisher_test <- mafCompare(m1 = data_quiescent_temporary_maf, m2 = data_fast.cycling_temporary_maf, 
                              m1Name = paste("HQ,", as.character(data_merged_pancancer$CancerType)[1]), 
                              m2Name = paste("FC,", as.character(data_merged_pancancer$CancerType)[1]))
    
    fisher_test_results[[as.character(data_merged_pancancer$CancerType)[1]]] <- as.data.frame(fisher_test[[1]])
    
    
    if (any(fisher_test[[1]]$adjPval < 0.05) == TRUE){
      print(forestPlot(fisher_test, fdr= 0.05, titleSize = 1.6, geneFontSize = 1.1, lineWidth = 1.3))
    } else {
      print(c("No significant genes for this cancer", cancer_files[i]))
      cancers_no_sig_genes <- c(cancers_no_sig_genes, as.character(data_merged_pancancer$CancerType)[1])
    }
      
    } else {
      print(c("Not enough Quiescent or Fast Cycling samples for this cancer", cancer_files[i]))
      cancers_not_enough_samples <- c(cancers_not_enough_samples, as.character(data_merged_pancancer$CancerType)[1])
    }
    
    } else {
      print(c("No Barcode match for this cancer", cancer_files[i]))
      cancers_no_barcode_match <- c(cancers_no_barcode_match, as.character(data_merged_pancancer$CancerType)[1])
    }
}

#making plots with sample size for each cancer
sample_size_dataframe <- na.omit(sample_size_dataframe)

sample_size_dataframe_plot <- rbind(data.frame(Cancer = sample_size_dataframe[,1], Size = as.numeric(sample_size_dataframe[,2]), Group = "Highly Quiescent"),
                                    data.frame(Cancer = sample_size_dataframe[,1], Size = as.numeric(sample_size_dataframe[,3]), Group = "Fast Cycling"))


setwd("C:/Users/victo/Desktop/Fourth year project/Data/Maftools")
pdf(file= "balance_maftools_pancancer.pdf")
for (i in unique(sample_size_dataframe_plot$Cancer)){
  print(ggplot(data= sample_size_dataframe_plot[sample_size_dataframe_plot$Cancer == i,], aes(x="", y= Size, fill= Group))+ 
    geom_bar(stat="identity")+
    geom_col(color = "black") +
    geom_text(aes(label = Size), position = position_stack(vjust = 0.5))+
    coord_polar("y")+
    theme(axis.text = element_blank(), plot.title = element_text(hjust = 0.5))+
    theme_void()+
    ggtitle(i))
}
dev.off()

#BRCA
#reading the BRCA MAF file
setwd("C:/Users/victo/Desktop/Fourth year project/Data/Maftools/GDCdata_SNV")
brca <- read.maf("BRCA.gz")
brca <- data.frame(brca@data)

#converting the Barcodes for BRCA TCGA to characters and splitting the Barcode
brca$Barcode2 <- sapply(as.character(brca$Tumor_Sample_Barcode), function(x) paste(strsplit(x,"-")[[1]][1:3],collapse="-"))

#reading clinical data for BRCA
brca_clinical <- GDCquery_clinic("TCGA-BRCA", type= "clinical", save.csv= FALSE)
#merging the clinical data with pancancer
brca_clinical <- merge(brca_clinical, pancancer_quiescence_groups, by.x= "submitter_id", by.y= "Barcode2")

#merging the brca_data and pancancer dataframes according to Barcode 2 
brca_merged_pancancer <- merge(brca, pancancer_quiescence_groups, by= "Barcode2")

#separating the data into two dataframes: for Highly Quiescent + Fast Cycling
brca_quiescent <- brca_merged_pancancer[brca_merged_pancancer$Group == "Highly Quiescent",]
brca_fast.cycling <- brca_merged_pancancer[brca_merged_pancancer$Group == "Fast Cycling",] 
#reading the two dataframes as maf files
brca_quiescent.maf <- read.maf(brca_quiescent)
brca_fast.cycling.maf <- read.maf(brca_fast.cycling)
  
#oncoplots for quiescent and fast cycling
oncoplot(brca_quiescent.maf)
oncoplot(brca_fast.cycling.maf)

plotmafSummary(brca_quiescent.maf)
plotmafSummary(brca_fast.cycling.maf) 

OncogenicPathways(brca_quiescent.maf)
OncogenicPathways(brca_fast.cycling.maf)

#CoBar plot for brca cells
coBarplot(m1= brca_quiescent.maf, m2= brca_fast.cycling.maf, m1Name= "BRCA, Quiescent", m2Name= "BRCA, Fast Cycling")

#Fisher test to compare differentially mutated genes between quiescent and fast-cycling BRCA cells
brca_fisher_test <- mafCompare(m1= brca_quiescent.maf, m2= brca_fast.cycling.maf, m1Name= "BRCA, Quiescent", m2Name= "BRCA, Fast-cycling", minMut = 5)
forestPlot(mafCompareRes = brca_fisher_test, fdr = 0.05)
brca_fisher_test <- as.data.frame(brca_fisher_test[[1]])
  
plotOncodrive(oncodrive(brca_fast.cycling.maf))
plotOncodrive(oncodrive(brca_quiescent.maf))
  
somaticInteractions(maf = brca_quiescent.maf, pvalue = c(0.05, 0.1), fontSize = 0.5)
somaticInteractions(maf = brca_fast.cycling.maf, top = 25, pvalue = c(0.05, 0.1), fontSize = 0.5)

  
        #PI3K pathway analysis for BRCA
#checking genes in the PI3K pathway in BRCA quiescent and fast-cycling cells
PlotOncogenicPathways(maf = brca_quiescent.maf, pathways = "PI3K")
PlotOncogenicPathways(maf = brca_fast.cycling.maf, pathways = "PI3K")
  
#save the genes (added manually) from the pathway plots
PI3K_genes_quiescent <- c("PIK3CA", "AKT1", "PTEN", "PIK3R1", "MTOR", "RICTOR", "RPTOR", "AKT3", "DEPDC5", "PPP2R1A", "TSC2", "DEPTOR", "INPP4B", "MAPKAP1", "PIK3CB", "PIK3R2", "RPS6", "STK11")
PI3K_genes_fast.cycling <- c("PIK3CA", "PTEN", "PIK3R1", "MTOR", "AKT1", "TSC1", "PIK3CB", "INPP4B", "RPTOR", "DEPDC5", "RICTOR", "AKT3", "MAPKAP1", "PIK3R3", "RPS6", "RPS6KB1", "TSC2", "AKT2", "NPRL2", "PDK1", "PIK3R2", "PPP2R1A", "STK11", "DEPTOR", "NPRL3", "RHEB")

#Isolating samples that have at least one gene mutated in the PI3K pathway
#subsetting maf files for PI3K pathway genes
#for quiescent


#making an oncoplot for quiescent
oncoplot(brca_quiescent_PI3K.genes.selected_maf, top = as.numeric(length(PI3K_genes_quiescent)))
  
#for fast-cycling
brca_fast.cycling_PI3K.genes.selected_data <- data.frame()
for (i in 1:nrow(brca_fast.cycling)){
  if (brca_fast.cycling$Hugo_Symbol[i] %in% PI3K_genes_fast.cycling){
    brca_fast.cycling_PI3K.genes.selected_data <- rbind(brca_fast.cycling_PI3K.genes.selected_data, brca_data_fast.cycling[i,])
  }
}
brca_fast.cycling_PI3K.genes.selected_maf <- read.maf(maf= brca_fast.cycling_PI3K.genes.selected_data)
#making an oncoplot for fast-cycling
oncoplot(brca_fast.cycling_PI3K.genes.selected_maf, top = as.numeric(length(PI3K_genes_fast.cycling)))

#co oncoplot for quiescent and fast cycling cells with at least one mutation in one of the PI3K pathway gemes
coOncoplot(m1= brca_fast.cycling_PI3K.genes.selected_maf, m2= brca_quiescent_PI3K.genes.selected_maf, m1Name = "Fast-cycling, PI3K genes, BRCA", m2Name = "Quiescent, PI3K genes, BRCA", genes = PI3K_genes_fast.cycling)

#all quiescent genes found in fast-cycling
#print the genes found in fast-cycling but not quiescent 
print(PI3K_genes_fast.cycling[which(!PI3K_genes_fast.cycling %in% PI3K_genes_quiescent)])
  
#Fisher test for PI3K genes (between quiescent and fast-cycling), for the samples that have at least one mutation in one PI3K pathway gene
quiescent_vs_fast.cycling_brca_PI3K_fisher <- mafCompare(m1= brca_quiescent_PI3K.genes.selected_maf, m2= brca_fast.cycling_PI3K.genes.selected_maf, m1Name= "BRCA, PI3K genes, Quiescent", m2Name= "BRCA, PI3K genes, Fast-cycling", minMut = 5)
print(quiescent_vs_fast.cycling_brca_PI3K_fisher)
forestPlot(mafCompareRes = quiescent_vs_fast.cycling_brca_PI3K_fisher, fdr = 0.05)

#isolating Fisher test results with PI3K genes, considering all samples
fisher_results <- as.data.frame(quiescent_vs_fast.cycling_brca_fisher[[1]])
fisher_PI3K_genes.results <- data.frame()
for (i in 1:nrow(fisher_results)){
  if (fisher_results$Hugo_Symbol[i] %in% PI3K_genes_quiescent | fisher_results$Hugo_Symbol[i] %in% PI3K_genes_fast.cycling){
    fisher_PI3K_genes.results <- rbind(fisher_PI3K_genes.results, fisher_results[i,])
  }
}
print(fisher_PI3K_genes.results)

    #comparing quiescence score with stage of cancer, BRCA
#fast-cycling and quiescent
res_t.test_quiescence_stage_brca <- brca_clinical %>% t_test(QuiescenceScore  ~ ajcc_pathologic_stage)
write.csv2(res_t.test_quiescence_stage_brca, "/Users/victorkonstantellos/Desktop/Fourth year project/Data/Maftools/res_t_test_quiescence_stage_brca.csv", col.names =  TRUE)
brca_all.data_merged_pancancer$ajcc_pathologic_stage <- factor(brca_all.data_merged_pancancer$ajcc_pathologic_stage, levels = c("Stage I", "Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB", "Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV", "Stage X", NA))
ggboxplot(data = brca_all.data_merged_pancancer, x="ajcc_pathologic_stage", y="QuiescenceScore", title = "Quiescence ~ Stage, BRCA", color= "Group")
res_anova_quiescence_stage_brca<- aov(QuiescenceScore ~ ajcc_pathologic_stage, data = brca_all.data_merged_pancancer)
summary.aov(res_anova_quiescence_stage_brca)

#merging stages numbers together
brca_clinical$ajcc_pathologic_stage_number <- rep(NA, nrow(brca_clinical))
for (i in 1:nrow(brca_clinical)){
  if (is.na(brca_clinical$ajcc_pathologic_stage[i]) == FALSE & !brca_clinical$ajcc_pathologic_stage[i] == "Stage X"){
    if (brca_clinical$ajcc_pathologic_stage[i] == "Stage I" || brca_clinical$ajcc_pathologic_stage[i] == "Stage IA" || brca_clinical$ajcc_pathologic_stage[i] == "Stage IB"){
      brca_clinical$ajcc_pathologic_stage_number[i] <- "I"
    } 
    if (brca_clinical$ajcc_pathologic_stage[i] == "Stage II" || brca_clinical$ajcc_pathologic_stage[i] == "Stage IIA" || brca_clinical$ajcc_pathologic_stage[i] == "Stage IIB"){
      brca_clinical$ajcc_pathologic_stage_number[i] <- "II"
    }
    if (brca_clinical$ajcc_pathologic_stage[i] == "Stage III" || brca_clinical$ajcc_pathologic_stage[i] == "Stage IIIA" || brca_clinical$ajcc_pathologic_stage[i] == "Stage IIIB" || brca_clinical$ajcc_pathologic_stage[i] == "Stage IIIC"){
      brca_clinical$ajcc_pathologic_stage_number[i] <- "III"
    }
    if (brca_clinical$ajcc_pathologic_stage[i] == "Stage IV"){
      brca_clinical$ajcc_pathologic_stage_number[i] <- "IV"
    } 
  }
}
res_t.test_quiescence_stage_number_brca <- brca_clinical %>% t_test(QuiescenceScore  ~ ajcc_pathologic_stage_number)
brca_clinical$ajcc_pathologic_stage_number <- factor(brca_clinical$ajcc_pathologic_stage_number, levels = c("I", "II", "III", "IV", NA))
ggboxplot(data = brca_clinical, x="ajcc_pathologic_stage_number", y="QuiescenceScore", title = "Quiescence ~ Stage, BRCA", color= "Group")

    #comparing quiescence score with age of diagnosis, fast-cycling and quiescent, BRCA
ggscatter(data= brca_clinical, x="age_at_diagnosis", y="QuiescenceScore", title = "Quiescence ~ Age at Diagnosis, BRCA", color = "Group", add = "reg.line")+ stat_cor(aes(color = Group), label.x = 11000)

    #comparing quiescence score with prior malignancy
#fast-cycling and quiescent
res_t.test_quiescence_prior.malignancy_brca <- brca_clinical %>% t_test(QuiescenceScore  ~ prior_malignancy)
ggboxplot(data = brca_clinical, x="prior_malignancy", y="QuiescenceScore", title = "Quiescence ~ Prior Malignancy, BRCA")+ stat_pvalue_manual(res_t.test_quiescence_prior.malignancy_brca, hide.ns= TRUE, y.position = c(13,15,17))

    #comparing quiescence score with race
res_t.test_quiescence_race_brca <- brca_clinical %>% t_test(QuiescenceScore  ~ race)
ggboxplot(data = brca_clinical, x="race", y="QuiescenceScore", title = "Quiescence ~ Race, BRCA", color= "Group")+ stat_pvalue_manual(res_t.test_quiescence_race_brca, hide.ns= TRUE, y.position = c(10,14))

#checking which cancer types are found in pancancer
if (askYesNo("Do you want to find which of the cancer types match with pancancer?") == TRUE){
  found_in_pancancer <- c()
  for (i in cancer_types){
    maf <- GDCquery_Maf(i, pipelines = "mutect2") %>% read.maf
    df <- data.frame(maf@data)
    df$Tumor_Sample_Barcode <- as.character(df$Tumor_Sample_Barcode)
    for (j in 1:nrow(df)){
      a <- strsplit(df$Tumor_Sample_Barcode[j], "-")
      df$Barcode2[j]= paste(a[[1]][1:3], collapse = "-")
    }
    if(any(df$Barcode2 %in% pancancer_quiescence_groups$Barcode2) == TRUE){
      found_in_pancancer <- c(found_in_pancancer, i)
    }
  }
}

#"DLBC" and "LAML" not in pancancer

#finding which of the cancer types that are found in pancancer can be plotted with oncoplot
if(askYesNo("Do you want to determine which cancer types that are found in the pancancer work with maftools oncoplots?") == TRUE){
  found_in_pancancer_and_work <- c()
  for (i in found_in_pancancer){
    maf <- GDCquery_Maf(i, pipelines = "mutect2") %>% read.maf
    df <- data.frame(maf@data)
    df$CancerStudy <- i
    df$Tumor_Sample_Barcode <- as.character(df$Tumor_Sample_Barcode)
    for (j in 1:nrow(df)){
      a <- strsplit(df$Tumor_Sample_Barcode[j], "-")
      df$Barcode2[j]= paste(a[[1]][1:3], collapse = "-")
    }
    df_merged_pancancer <- merge(df, pancancer_quiescence_groups, by.x = "Barcode2", by.y= "Barcode2")
    
    df_quiescent <- df_merged_pancancer[df_merged_pancancer$Group == "Highly Quiescent",]
    df_fast.cycling <- df_merged_pancancer[df_merged_pancancer$Group == "Fast Cycling",] 
    
    if (nrow(df_quiescent) & nrow(df_fast.cycling) > 0){
      found_in_pancancer_and_work <- c(found_in_pancancer_and_work, i)
    }
  }
  found_in_pancancer_and_work <- c("BLCA", "LGG",  "BRCA", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "UCEC", "UVM") 
  
}

#plotting
cancer_types <- found_in_pancancer_and_work
cancer_types_dataframe <- data.frame()
for (i in cancer_types){
  maf <- GDCquery_Maf(i, pipelines = "mutect2") %>% read.maf
  df <- data.frame(maf@data)
  df$CancerStudy <- i
  df$Tumor_Sample_Barcode <- as.character(df$Tumor_Sample_Barcode)
  for (j in 1:nrow(df)){
    a <- strsplit(df$Tumor_Sample_Barcode[j], "-")
    df$Barcode2[j]= paste(a[[1]][1:3], collapse = "-")
  }
  df_merged_pancancer <- merge(df, pancancer_quiescence_groups, by.x = "Barcode2", by.y= "Barcode2")
  
  df_quiescent <- df_merged_pancancer[df_merged_pancancer$Group == "Highly Quiescent",]
  df_fast.cycling <- df_merged_pancancer[df_merged_pancancer$Group == "Fast Cycling",] 
  
  df_quiescent.maf <- read.maf(df_quiescent)
  df_fast.cycling.maf <- read.maf(df_fast.cycling)
  
  oncoplot(df_quiescent.maf, top = 10)
  oncoplot(df_fast.cycling.maf, top = 10)
}

#reading the cancer_types_merged_pancancer
if(askYesNo("Do you want to read the cancer types merged with pancancer file?") == TRUE){
  cancer_types_merged_pancancer <- read.csv2("TCGA_data_merged_pancancer.csv")
}


#path to TCGA LAML MAF file
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
#clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 
laml = read.maf(maf = laml.maf, clinicalData = laml.clin)

brca <- GDCquery_Maf("BRCA", pipelines = "mutect2") %>% read.maf

brca.clinical <- GDCquery_clinic("TCGA-BRCA", type= "clinical", save.csv= FALSE)
for (i in 1:nrow(brca_data)){
  a <- strsplit(brca_data$Tumor_Sample_Barcode[i], "-")
  brca_data$Barcode2[i]= paste(a[[1]][1:3], collapse = "-")
}








