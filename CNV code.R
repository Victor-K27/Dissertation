#load the data
setwd("~/Desktop/Fourth year project/Data/CNV/segmentation files")
load("myEnvironment.RData")

#loading libraries
library(copynumber)
library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(rstatix)
library(cluster)
library(jaccard)
library(ggfortify)
library(mygene)
library(enrichR)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

hg38_genome <- read.delim("geneCoordinates_hg38.txt")
hg38_genome[hg38_genome$Chromosome.scaffold.name == "X", ]$Chromosome.scaffold.name <- 23
hg38_genome$segment <-paste(hg38_genome$Gene.start..bp., "-", hg38_genome$Gene.end..bp.)
#reading pancancer file and splitting barcode
pancancer_quiescence_groups <- read.csv2("PancancerQuiescenceGroupsCorrected.csv")
for (i in 1:nrow(pancancer_quiescence_groups)){
  a <- strsplit(pancancer_quiescence_groups$Barcode[i], "-")
  pancancer_quiescence_groups$Barcode2[i]= paste(a[[1]][1:4], collapse = "-")
}
#plotting the balancer of each group in all cancers in the pancancerquiescencegroupsfiles
#the data for the table was creating manually in excel
balance_all_cancers <- read.csv2("balance_all_cancers.csv")
ggplot(balance_all_cancers, aes(x= Cancer, y= Count,fill= Group))+ 
  geom_bar(stat = "identity")+
  ggtitle("Balance of Fast Cycling and Highly Quiescent cells for 31 cancers")
#plotting pie chart for sum
sum_balance_all_cancers <- data.frame(Count= c(sum(balance_all_cancers[balance_all_cancers$Group == "Highly Quiescent",]$Count),
                                                   sum(balance_all_cancers[balance_all_cancers$Group == "Fast Cycling",]$Count)), 
                                      Group= c("Fast Cycling", "Highly Quiescent"))
ggplot(sum_balance_all_cancers, aes(x= "", y = Count, fill = Group))+
  geom_col(color= "black") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5))+
  coord_polar(theta = "y")+
  ggtitle("Total samples in 31 cancers")
#plotting pie charts for all 31 cancers individually
for (i in unique(pancancer_quiescence_groups$CancerType)){
  sum_dataframe_temp <- data.frame(Count= c(balance_all_cancers[balance_all_cancers$Cancer == i & balance_all_cancers$Group == "Highly Quiescent",]$Count,
                                            balance_all_cancers[balance_all_cancers$Cancer == i & balance_all_cancers$Group == "Fast Cycling",]$Count),
                                   Group= c("Fast Cycling", "Highly Quiescent"))
  print(ggplot(sum_dataframe_temp, aes(x= "", y = Count, fill = Group))+
    geom_col(color= "black") +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5))+
    coord_polar(theta = "y")+
    ggtitle(i))
}

#function to compare quiescent and fast cycling CNV in the same cancer
segmented_analysis <- function(cancer_type_segmentation_file){
  start_time <- Sys.time()
  cancer_type_segmentation_file <- cancer_type_segmentation_file[,-5]
  colnames(cancer_type_segmentation_file) <- c("sampleID", "chrom", "start.pos", "end.pos", "mean")
  #adding "A" at the end of barcode for ovarian data, this was missing
  cancer_type_segmentation_file$sampleID <- gsub(" ", "", paste(cancer_type_segmentation_file$sampleID, "A"))
  
  cn_merged_pancancer <- merge(x=cancer_type_segmentation_file, y=pancancer_quiescence_groups, by.x= "sampleID", by.y= "Barcode2")
  
  if (length(unique(cn_merged_pancancer$Group)) > 1){
    cn_quiescent <- cn_merged_pancancer[which(cn_merged_pancancer$Group == "Highly Quiescent"),]
    cn_fast.cycling <- cn_merged_pancancer[which(cn_merged_pancancer$Group == "Fast Cycling"),]
    
    test_dataframe <- data.frame()
    test_dataframe_significant <- data.frame()
    overlap_dataframe_aggregated_all <- data.frame()
    for (i in unique(cancer_type_segmentation_file$chrom)[1]){
      seq <- seq(from = min(hg38_genome[hg38_genome$Chromosome.scaffold.name == i,]$Gene.start..bp.), to = max(hg38_genome[hg38_genome$Chromosome.scaffold.name == i,]$Gene.end..bp.), by = 100000)
      sequence <- c()
      for (j in 1:length(seq)-1){
        sequence <- c(sequence, paste(seq[j], "-", seq[j+1]))
      }
      sequence <- sequence[-1]
      
      gr_hg38_genome <- GRanges(
        seqnames = Rle(i),
        ranges = IRanges(sequence))
      
      cn_quiescent_samples <- cn_quiescent[cn_quiescent$chrom == i,]
      gr_quiescent <- GRanges(
        seqnames = Rle(cn_quiescent_samples$chrom),
        ranges = IRanges(cn_quiescent_samples$start.pos, end= cn_quiescent_samples$end.pos, names = cn_quiescent_samples$sampleID),
        Mean = cn_quiescent_samples$mean,
        Group = cn_quiescent_samples$Group,
        SampleID = cn_quiescent_samples$sampleID)
      
      cn_fast.cycling_samples <- cn_fast.cycling[cn_fast.cycling$chrom == i,]
      gr_fast.cycling <- GRanges(
        seqnames = Rle(cn_fast.cycling_samples$chrom),
        ranges = IRanges(cn_fast.cycling_samples$start.pos, end= cn_fast.cycling_samples$end.pos, names = cn_fast.cycling_samples$sampleID),
        Mean = cn_fast.cycling_samples$mean,
        Group = cn_fast.cycling_samples$Group,
        SampleID = cn_fast.cycling_samples$sampleID)
      
      overlaps_quiescent <- as.data.frame(findOverlaps(gr_hg38_genome, gr_quiescent))
      overlap_dataframe_quiescent <- cbind(as.data.frame(gr_hg38_genome[overlaps_quiescent$queryHits], row.names = seq(from = 1, to = nrow(overlaps_quiescent), 1)), as.data.frame(gr_quiescent[overlaps_quiescent$subjectHits], row.names = seq(from = 1, to = nrow(overlaps_quiescent), by = 1)))
      colnames(overlap_dataframe_quiescent) <- c("chr_1", "start_1", "end_1", "width_1", "strand_1", "seq_2", "start_2", "end_2", "width_2", "strand_2", "mean", "group", "sampleID")
      
      overlaps_fast.cycling <- as.data.frame(findOverlaps(gr_hg38_genome, gr_fast.cycling))
      overlap_dataframe_fast.cycling <- cbind(as.data.frame(gr_hg38_genome[overlaps_fast.cycling$queryHits], row.names = seq(from = 1, to = nrow(overlaps_fast.cycling), by = 1)), as.data.frame(gr_fast.cycling[overlaps_fast.cycling$subjectHits], row.names = seq(from = 1, to = nrow(overlaps_fast.cycling), by = 1)))
      colnames(overlap_dataframe_fast.cycling) <- c("chr_1", "start_1", "end_1", "width_1", "strand_1", "seq_2", "start_2", "end_2", "width_2", "strand_2", "mean", "group", "sampleID")
      
      overlap_dataframe_quiescent$segment <- paste(overlap_dataframe_quiescent$start_1, "-", overlap_dataframe_quiescent$end_1)
      quiecent_aggregated <- aggregate(data = overlap_dataframe_quiescent, mean ~ segment, mean)
      quiecent_aggregated$count <- aggregate(data = overlap_dataframe_quiescent, mean ~ segment, FUN = length)[,2]
      quiecent_aggregated$sd <- aggregate(data = overlap_dataframe_quiescent, mean ~ segment, sd)[,2]
      quiecent_aggregated$mean_plus_sd <- as.numeric(quiecent_aggregated$mean)+as.numeric(quiecent_aggregated$sd)
      quiecent_aggregated$mean_minus_sd <- as.numeric(quiecent_aggregated$mean)-as.numeric(quiecent_aggregated$sd)
      colnames(quiecent_aggregated) <- c("Segment", "Mean", "Count", "SD", "Mean+SD", "Mean-SD")
      quiecent_aggregated$group <- "Highly Quiescent"
      quiecent_aggregated$cancer_type <- cn_merged_pancancer$CancerType[1]
      quiecent_aggregated$start.pos <- as.numeric(vapply(strsplit(quiecent_aggregated$Segment,"-"), `[`, 1, FUN.VALUE=character(1)))
      quiecent_aggregated$end.pos <- as.numeric(vapply(strsplit(quiecent_aggregated$Segment,"-"), `[`, 2, FUN.VALUE=character(1)))
      
      overlap_dataframe_fast.cycling$segment <- paste(overlap_dataframe_fast.cycling$start_1, "-", overlap_dataframe_fast.cycling$end_1)
      fast.cycling_aggregated <- aggregate(data = overlap_dataframe_fast.cycling, mean ~ segment, mean)
      fast.cycling_aggregated$count <- aggregate(data = overlap_dataframe_fast.cycling, mean ~ segment, FUN = length)[,2]
      fast.cycling_aggregated$sd <- aggregate(data = overlap_dataframe_fast.cycling, mean ~ segment, sd)[,2]
      fast.cycling_aggregated$mean_plus_sd <- as.numeric(fast.cycling_aggregated$mean)+as.numeric(fast.cycling_aggregated$sd)
      fast.cycling_aggregated$mean_minus_sd <- as.numeric(fast.cycling_aggregated$mean)-as.numeric(fast.cycling_aggregated$sd)
      colnames(fast.cycling_aggregated) <- c("Segment", "Mean", "Count", "SD", "Mean+SD", "Mean-SD")
      fast.cycling_aggregated$group <- "Fast Cycling"
      fast.cycling_aggregated$cancer_type <- cn_merged_pancancer$CancerType[1]
      fast.cycling_aggregated$start.pos <- as.numeric(vapply(strsplit(fast.cycling_aggregated$Segment,"-"), `[`, 1, FUN.VALUE=character(1)))
      fast.cycling_aggregated$end.pos <- as.numeric(vapply(strsplit(fast.cycling_aggregated$Segment,"-"), `[`, 2, FUN.VALUE=character(1)))
      
      overlap_dataframe_aggregated <- rbind(quiecent_aggregated, fast.cycling_aggregated)
      overlap_dataframe_aggregated$chrom <- i
      overlap_dataframe_aggregated_all <- rbind(overlap_dataframe_aggregated_all, overlap_dataframe_aggregated)
      overlap_dataframe <- rbind(overlap_dataframe_quiescent, overlap_dataframe_fast.cycling)
      
      test <- overlap_dataframe %>% group_by(segment) %>% wilcox_test(mean ~ group) %>% adjust_pvalue(method = "bonferroni") %>% add_significance("p.adj")
      test$chrom <- i
      test$cancer_type <- cn_merged_pancancer$CancerType[1]
      test_dataframe <- rbind(test_dataframe, test)
    }
    test_dataframe_significant <- rbind(test_dataframe_significant, test_dataframe[test_dataframe$p.adj <= 0.05,])
    return(list(test_dataframe, test_dataframe_significant, overlap_dataframe_aggregated_all))
  } else {
    return(paste("No HQ or FC for this cancer", cn_merged_pancancer$CancerType[1]))
  }
  print(Sys.time()- start_time)
}
#function to compare HQ ~ HQ, FC ~ FC, HQ ~ FC, and FC ~ HQ between two cancers 
segmented_analysis_2 <- function(cancer_type_segmentation_file_1, cancer_type_segmentation_file_2, window_width, comparisons){
  start_time <- Sys.time()
  cn_merged_pancancer_1 <- merge(x=cancer_type_segmentation_file_1, y=pancancer_quiescence_groups, by.x= "sampleID", by.y= "Barcode2")
  cn_merged_pancancer_2 <- merge(x=cancer_type_segmentation_file_2, y=pancancer_quiescence_groups, by.x= "sampleID", by.y= "Barcode2")
  
  cn_quiescent_1 <- cn_merged_pancancer_1[cn_merged_pancancer_1$Group == "Highly Quiescent",]
  cn_fast.cycling_1 <- cn_merged_pancancer_1[cn_merged_pancancer_1$Group == "Fast Cycling",]
  cn_quiescent_2 <- cn_merged_pancancer_2[cn_merged_pancancer_2$Group == "Highly Quiescent",]
  cn_fast.cycling_2 <- cn_merged_pancancer_2[cn_merged_pancancer_2$Group == "Fast Cycling",]
  
  test_dataframe <- data.frame()
  overlap_dataframe_aggregated_all <- data.frame()
  chromosome_selection <- 
    for (i in unique(cn_merged_pancancer_2$chrom)){
      seq <- seq(from = min(hg38_genome[hg38_genome$Chromosome.scaffold.name == i,]$Gene.start..bp.), to = max(hg38_genome[hg38_genome$Chromosome.scaffold.name == i,]$Gene.end..bp.), by = window_width)
      sequence <- c()
      for (j in 1:length(seq)-1){
        sequence <- c(sequence, paste(seq[j], "-", seq[j+1]))
      }
      sequence <- sequence[-1]
      
      gr_hg38_genome <- GRanges(
        seqnames = Rle(i),
        ranges = IRanges(sequence))
      
      cn_quiescent_samples_1 <- cn_quiescent_1[cn_quiescent_1$chrom == i,]
      gr_quiescent_1 <- GRanges(
        seqnames = Rle(cn_quiescent_samples_1$chrom),
        ranges = IRanges(cn_quiescent_samples_1$start.pos, end= cn_quiescent_samples_1$end.pos, names = cn_quiescent_samples_1$sampleID),
        Mean = cn_quiescent_samples_1$mean,
        Group = cn_quiescent_samples_1$Group,
        SampleID = cn_quiescent_samples_1$sampleID)
      
      cn_quiescent_samples_2 <- cn_quiescent_2[cn_quiescent_2$chrom == i,]
      gr_quiescent_2 <- GRanges(
        seqnames = Rle(cn_quiescent_samples_2$chrom),
        ranges = IRanges(cn_quiescent_samples_2$start.pos, end= cn_quiescent_samples_2$end.pos, names = cn_quiescent_samples_2$sampleID),
        Mean = cn_quiescent_samples_2$mean,
        Group = cn_quiescent_samples_2$Group,
        SampleID = cn_quiescent_samples_2$sampleID)
      
      cn_fast.cycling_samples_1 <- cn_fast.cycling_1[cn_fast.cycling_1$chrom == i,]
      gr_fast.cycling_1 <- GRanges(
        seqnames = Rle(cn_fast.cycling_samples_1$chrom),
        ranges = IRanges(cn_fast.cycling_samples_1$start.pos, end= cn_fast.cycling_samples_1$end.pos, names = cn_fast.cycling_samples_1$sampleID),
        Mean = cn_fast.cycling_samples_1$mean,
        Group = cn_fast.cycling_samples_1$Group,
        SampleID = cn_fast.cycling_samples_1$sampleID)
      
      cn_fast.cycling_samples_2 <- cn_fast.cycling_2[cn_fast.cycling_2$chrom == i,]
      gr_fast.cycling_2 <- GRanges(
        seqnames = Rle(cn_fast.cycling_samples_2$chrom),
        ranges = IRanges(cn_fast.cycling_samples_2$start.pos, end= cn_fast.cycling_samples_2$end.pos, names = cn_fast.cycling_samples_2$sampleID),
        Mean = cn_fast.cycling_samples_2$mean,
        Group = cn_fast.cycling_samples_2$Group,
        SampleID = cn_fast.cycling_samples_2$sampleID)
      
      overlaps_quiescent_1 <- as.data.frame(findOverlaps(gr_hg38_genome, gr_quiescent_1))
      overlap_dataframe_quiescent_1 <- cbind(as.data.frame(gr_hg38_genome[overlaps_quiescent_1$queryHits], row.names = seq(from = 1, to = nrow(overlaps_quiescent_1), 1)), as.data.frame(gr_quiescent_1[overlaps_quiescent_1$subjectHits], row.names = seq(from = 1, to = nrow(overlaps_quiescent_1), by = 1)))
      colnames(overlap_dataframe_quiescent_1) <- c("chr_1", "start_1", "end_1", "width_1", "strand_1", "seq_2", "start_2", "end_2", "width_2", "strand_2", "mean", "group", "sampleID")
      
      overlaps_quiescent_2 <- as.data.frame(findOverlaps(gr_hg38_genome, gr_quiescent_2))
      overlap_dataframe_quiescent_2 <- cbind(as.data.frame(gr_hg38_genome[overlaps_quiescent_2$queryHits], row.names = seq(from = 1, to = nrow(overlaps_quiescent_2), 1)), as.data.frame(gr_quiescent_2[overlaps_quiescent_2$subjectHits], row.names = seq(from = 1, to = nrow(overlaps_quiescent_2), by = 1)))
      colnames(overlap_dataframe_quiescent_2) <- c("chr_1", "start_1", "end_1", "width_1", "strand_1", "seq_2", "start_2", "end_2", "width_2", "strand_2", "mean", "group", "sampleID")
      
      overlaps_fast.cycling_1 <- as.data.frame(findOverlaps(gr_hg38_genome, gr_fast.cycling_1))
      overlap_dataframe_fast.cycling_1 <- cbind(as.data.frame(gr_hg38_genome[overlaps_fast.cycling_1$queryHits], row.names = seq(from = 1, to = nrow(overlaps_fast.cycling_1), by = 1)), as.data.frame(gr_fast.cycling_1[overlaps_fast.cycling_1$subjectHits], row.names = seq(from = 1, to = nrow(overlaps_fast.cycling_1), by = 1)))
      colnames(overlap_dataframe_fast.cycling_1) <- c("chr_1", "start_1", "end_1", "width_1", "strand_1", "seq_2", "start_2", "end_2", "width_2", "strand_2", "mean", "group", "sampleID")
      
      overlaps_fast.cycling_2 <- as.data.frame(findOverlaps(gr_hg38_genome, gr_fast.cycling_2))
      overlap_dataframe_fast.cycling_2 <- cbind(as.data.frame(gr_hg38_genome[overlaps_fast.cycling_2$queryHits], row.names = seq(from = 1, to = nrow(overlaps_fast.cycling_2), by = 1)), as.data.frame(gr_fast.cycling_2[overlaps_fast.cycling_2$subjectHits], row.names = seq(from = 1, to = nrow(overlaps_fast.cycling_2), by = 1)))
      colnames(overlap_dataframe_fast.cycling_2) <- c("chr_1", "start_1", "end_1", "width_1", "strand_1", "seq_2", "start_2", "end_2", "width_2", "strand_2", "mean", "group", "sampleID")
      
      overlap_dataframe_quiescent_1$segment <- paste(overlap_dataframe_quiescent_1$start_1, "-", overlap_dataframe_quiescent_1$end_1)
      overlap_dataframe_quiescent_1$cancer_type <- cn_merged_pancancer_1$CancerType[1]
      quiecent_aggregated_1 <- aggregate(data = overlap_dataframe_quiescent_1, mean ~ segment, mean)
      quiecent_aggregated_1$count <- aggregate(data = overlap_dataframe_quiescent_1, mean ~ segment, FUN = length)[,2]
      quiecent_aggregated_1$sd <- aggregate(data = overlap_dataframe_quiescent_1, mean ~ segment, sd)[,2]
      quiecent_aggregated_1$mean_plus_sd <- as.numeric(quiecent_aggregated_1$mean)+as.numeric(quiecent_aggregated_1$sd)
      quiecent_aggregated_1$mean_minus_sd <- as.numeric(quiecent_aggregated_1$mean)-as.numeric(quiecent_aggregated_1$sd)
      colnames(quiecent_aggregated_1) <- c("Segment", "Mean", "Count", "SD", "Mean+SD", "Mean-SD")
      quiecent_aggregated_1$group <- "Highly Quiescent"
      quiecent_aggregated_1$start.pos <- as.numeric(vapply(strsplit(quiecent_aggregated_1$Segment,"-"), `[`, 1, FUN.VALUE=character(1)))
      quiecent_aggregated_1$end.pos <- as.numeric(vapply(strsplit(quiecent_aggregated_1$Segment,"-"), `[`, 2, FUN.VALUE=character(1)))
      
      overlap_dataframe_quiescent_2$segment <- paste(overlap_dataframe_quiescent_2$start_1, "-", overlap_dataframe_quiescent_2$end_1)
      overlap_dataframe_quiescent_2$cancer_type <- cn_merged_pancancer_2$CancerType[1]
      quiecent_aggregated_2 <- aggregate(data = overlap_dataframe_quiescent_2, mean ~ segment, mean)
      quiecent_aggregated_2$count <- aggregate(data = overlap_dataframe_quiescent_2, mean ~ segment, FUN = length)[,2]
      quiecent_aggregated_2$sd <- aggregate(data = overlap_dataframe_quiescent_2, mean ~ segment, sd)[,2]
      quiecent_aggregated_2$mean_plus_sd <- as.numeric(quiecent_aggregated_2$mean)+as.numeric(quiecent_aggregated_2$sd)
      quiecent_aggregated_2$mean_minus_sd <- as.numeric(quiecent_aggregated_2$mean)-as.numeric(quiecent_aggregated_2$sd)
      colnames(quiecent_aggregated_2) <- c("Segment", "Mean", "Count", "SD", "Mean+SD", "Mean-SD")
      quiecent_aggregated_2$group <- "Highly Quiescent"
      quiecent_aggregated_2$start.pos <- as.numeric(vapply(strsplit(quiecent_aggregated_2$Segment,"-"), `[`, 1, FUN.VALUE=character(1)))
      quiecent_aggregated_2$end.pos <- as.numeric(vapply(strsplit(quiecent_aggregated_2$Segment,"-"), `[`, 2, FUN.VALUE=character(1)))
      
      overlap_dataframe_fast.cycling_1$segment <- paste(overlap_dataframe_fast.cycling_1$start_1, "-", overlap_dataframe_fast.cycling_1$end_1)
      overlap_dataframe_fast.cycling_1$cancer_type <- cn_merged_pancancer_1$CancerType[1]
      fast.cycling_aggregated_1 <- aggregate(data = overlap_dataframe_fast.cycling_1, mean ~ segment, mean)
      fast.cycling_aggregated_1$count <- aggregate(data = overlap_dataframe_fast.cycling_1, mean ~ segment, FUN = length)[,2]
      fast.cycling_aggregated_1$sd <- aggregate(data = overlap_dataframe_fast.cycling_1, mean ~ segment, sd)[,2]
      fast.cycling_aggregated_1$mean_plus_sd <- as.numeric(fast.cycling_aggregated_1$mean)+as.numeric(fast.cycling_aggregated_1$sd)
      fast.cycling_aggregated_1$mean_minus_sd <- as.numeric(fast.cycling_aggregated_1$mean)-as.numeric(fast.cycling_aggregated_1$sd)
      colnames(fast.cycling_aggregated_1) <- c("Segment", "Mean", "Count", "SD", "Mean+SD", "Mean-SD")
      fast.cycling_aggregated_1$group <- "Fast Cycling"
      fast.cycling_aggregated_1$start.pos <- as.numeric(vapply(strsplit(fast.cycling_aggregated_1$Segment,"-"), `[`, 1, FUN.VALUE=character(1)))
      fast.cycling_aggregated_1$end.pos <- as.numeric(vapply(strsplit(fast.cycling_aggregated_1$Segment,"-"), `[`, 2, FUN.VALUE=character(1)))
      
      overlap_dataframe_fast.cycling_2$segment <- paste(overlap_dataframe_fast.cycling_2$start_1, "-", overlap_dataframe_fast.cycling_2$end_1)
      overlap_dataframe_fast.cycling_2$cancer_type <- cn_merged_pancancer_2$CancerType[1]
      fast.cycling_aggregated_2 <- aggregate(data = overlap_dataframe_fast.cycling_2, mean ~ segment, mean)
      fast.cycling_aggregated_2$count <- aggregate(data = overlap_dataframe_fast.cycling_2, mean ~ segment, FUN = length)[,2]
      fast.cycling_aggregated_2$sd <- aggregate(data = overlap_dataframe_fast.cycling_2, mean ~ segment, sd)[,2]
      fast.cycling_aggregated_2$mean_plus_sd <- as.numeric(fast.cycling_aggregated_2$mean)+as.numeric(fast.cycling_aggregated_2$sd)
      fast.cycling_aggregated_2$mean_minus_sd <- as.numeric(fast.cycling_aggregated_2$mean)-as.numeric(fast.cycling_aggregated_2$sd)
      colnames(fast.cycling_aggregated_2) <- c("Segment", "Mean", "Count", "SD", "Mean+SD", "Mean-SD")
      fast.cycling_aggregated_2$group <- "Fast Cycling"
      fast.cycling_aggregated_2$start.pos <- as.numeric(vapply(strsplit(fast.cycling_aggregated_2$Segment,"-"), `[`, 1, FUN.VALUE=character(1)))
      fast.cycling_aggregated_2$end.pos <- as.numeric(vapply(strsplit(fast.cycling_aggregated_2$Segment,"-"), `[`, 2, FUN.VALUE=character(1)))
      
      if (comparisons == "HQ ~ HQ"){
        overlap_dataframe <- rbind(overlap_dataframe_quiescent_1[overlap_dataframe_quiescent_1$segment %in% overlap_dataframe_quiescent_2$segment,], overlap_dataframe_quiescent_2[overlap_dataframe_quiescent_2$segment %in% overlap_dataframe_quiescent_1$segment,])
        test <- overlap_dataframe %>% group_by(segment) %>% t_test(mean ~ cancer_type) %>% adjust_pvalue(method = "bonferroni") %>% add_significance("p.adj")
        test$chrom <- i
        test_dataframe <- rbind(test_dataframe, test)
        aggregated_dataframe <- rbind(quiecent_aggregated_1, quiecent_aggregated_2)
        overlap_dataframe_aggregated_all <- rbind(overlap_dataframe_aggregated_all, aggregated_dataframe)
      }
      if (comparisons == "FC ~ FC"){
        overlap_dataframe <- rbind(overlap_dataframe_fast.cycling_1[overlap_dataframe_fast.cycling_1$segment %in% overlap_dataframe_fast.cycling_2$segment,], overlap_dataframe_fast.cycling_2[overlap_dataframe_fast.cycling_2$segment %in% overlap_dataframe_fast.cycling_1$segment,])
        test <- overlap_dataframe %>% group_by(segment) %>% t_test(mean ~ cancer_type) %>% adjust_pvalue(method = "bonferroni") %>% add_significance("p.adj")
        test$chrom <- i
        test_dataframe <- rbind(test_dataframe, test)
        aggregated_dataframe <- rbind(fast.cycling_aggregated_1, fast.cycling_aggregated_2)
        overlap_dataframe_aggregated_all <- rbind(overlap_dataframe_aggregated_all, aggregated_dataframe)
      }
      if (comparisons == "HQ ~ FC"){
        overlap_dataframe <- rbind(overlap_dataframe_quiescent_1[overlap_dataframe_quiescent_1$segment %in% overlap_dataframe_fast.cycling_2$segment], overlap_dataframe_fast.cycling_2[overlap_dataframe_fast.cycling_2$segment %in% overlap_dataframe_quiescent_1$segment,])
        test <- overlap_dataframe %>% group_by(segment) %>% t_test(mean ~ group) %>% adjust_pvalue(method = "bonferroni") %>% add_significance("p.adj")
        test$chrom <- i
        test_dataframe <- rbind(test_dataframe, test)
        aggregated_dataframe <- rbind(quiecent_aggregated_1, fast.cycling_aggregated_2)
        overlap_dataframe_aggregated_all <- rbind(overlap_dataframe_aggregated_all, aggregated_dataframe)
      }
      if (comparisons == "FC ~ HQ"){
        overlap_dataframe <- rbind(overlap_dataframe_fast.cycling_1[overlap_dataframe_fast.cycling_1$segment %in% overlap_dataframe_quiescent_2$segment,], overlap_dataframe_quiescent_2[overlap_dataframe_quiescent_2$segment %in% overlap_dataframe_fast.cycling_1$segment,])
        test <- overlap_dataframe %>% group_by(segment) %>% t_test(mean ~ group) %>% adjust_pvalue(method = "bonferroni") %>% add_significance("p.adj")
        test$chrom <- i
        test_dataframe <- rbind(test_dataframe, test)
        aggregated_dataframe <- rbind(fast.cycling_aggregated_2, quiecent_aggregated_2)
        overlap_dataframe_aggregated_all <- rbind(overlap_dataframe_aggregated_all, aggregated_dataframe)
      }
    }
  return(list(test_dataframe, overlap_dataframe_aggregated_all))
  print(Sys.time()- start_time)
}

#running the analysis on 29 cancer types, takes about 7 hours
setwd("~/Desktop/Fourth year project/Data/CNV/segmentation files")
file.list <- dir(pattern = "seg")
files <- lapply(file.list, read.delim)
res <- sapply(files, segmented_analysis)



#removing the cancer types without sufficient samples
res_refined <- res[-c(1,8,26,28)]

#isolating the cancer types that have significant segments, irrespective of CNV threshold (-0.33, +0.33)
cancer_types_sig <- data.frame()
for (i in 1:length(res_refined)){
  dataframe_test_sig <- res_refined[[i]][[2]]
  dataframe_aggregated <- res_refined[[i]][[3]]
  dataframe_aggregated_sig <- dataframe_aggregated[dataframe_aggregated$Segment %in% dataframe_test_sig$segment,]
  test_sig_merged_aggregated_sig <- merge(dataframe_test_sig, dataframe_aggregated_sig, by.x= "segment", by.y= "Segment")
  test_sig_merged_aggregated_sig <- test_sig_merged_aggregated_sig[,-c(2:8,16,17,19,22)]
  
  cancer_types_sig <- rbind(cancer_types_sig, test_sig_merged_aggregated_sig)
}
  
#splitting by chromosomes
cancer_types_sig_chr <- cancer_types_sig %>% group_split(cancer_types_sig$chrom.x)

#creating the scoring dataframe
#changing the chr numbering (X -> 23)
hg38_genome$Chromosome.scaffold.name[hg38_genome$Chromosome.scaffold.name == "X"] <- 23
seg_chrom_all <- data.frame()
for (i in unique(hg38_genome$Chromosome.scaffold.name[hg38_genome$Chromosome.scaffold.name != "Y"])){
  seq <- seq(from = min(hg38_genome[hg38_genome$Chromosome.scaffold.name == i,]$Gene.start..bp.), to = max(hg38_genome[hg38_genome$Chromosome.scaffold.name == i,]$Gene.end..bp.), by = 100000)
  sequence <- c()
  for (j in 1:length(seq)-1){
    sequence <- c(sequence, paste(seq[j], "-", seq[j+1]))
  }
  sequence <- sequence[-1]
  seg_chrom <- data.frame(segment = sequence, chrom = i)
  seg_chrom_all <- rbind(seg_chrom_all, seg_chrom)
}
for (i in unique(cancer_types_sig$cancer_type.x)){
  seg_chrom_all[,i] <-  rep(0, nrow(seg_chrom_all))
}

#scoring the dataframe
for (i in unique(cancer_types_sig$chrom.x)){
  cancer_types_sig_temp <- cancer_types_sig[cancer_types_sig$chrom.x == i,]
  for (j in unique(cancer_types_sig_temp$cancer_type.x)){
    for (h in 1:nrow(seg_chrom_all[seg_chrom_all$chrom == i,])){
      if (seg_chrom_all[seg_chrom_all$chrom == i,]$segment[h] %in% cancer_types_sig_temp[cancer_types_sig_temp$cancer_type.x == j,]$segment){
        seg_chrom_all[seg_chrom_all$chrom == i,][h,j] <- 1
      } 
    }
  }
}
#adding the sum for each segment
seg_chrom_all["sum"] <- rep(0,nrow(seg_chrom_all))
for (i in 1:nrow(seg_chrom_all)){
  seg_chrom_all$sum[i] <- sum(seg_chrom_all[i,3:16])
}
#adding the sum for each cancer
colSums(seg_chrom_all[3:(ncol(seg_chrom_all)-1)])
sum_vector <- c()
for (i in 1:ncol(seg_chrom_all)){
  if (is.numeric(seg_chrom_all[,i]) == FALSE){
    sum_vector <- c(sum_vector, NA)
  } else {
    sum_vector <- c(sum_vector, sum(seg_chrom_all[,i]))
  }
}
seg_chrom_all <- rbind(seg_chrom_all, sum_vector)
#counting significant segments per chromosome
sum_chrom <- data.frame()
for (i in unique(seg_chrom_all$chrom)){
  chrom <- seg_chrom_all[seg_chrom_all$chrom == i,]
  data <- c(as.numeric(i), nrow(chrom))
  sum_chrom <- rbind(sum_chrom, data)
}
colnames(sum_chrom) <- c("Chrom", "sum")
sum_chrom <- sum_chrom[is.na(sum_chrom$Chrom) ==  FALSE,]

if(askYesNo("Do you want to currate the seg_chrom_all dataframe?") ==  TRUE){
  seg_chrom_all <- seg_chrom_all[,-c(6,7,11,16)]
}

#making new dataframes with the clusters
cluster_dataframe <- seg_chrom_all[,3:12]
cluster_dataframe <- cluster_dataframe[-nrow(cluster_dataframe),]
cluster_dataframe_2 <- as.data.frame(t(cluster_dataframe))
cluster_dataframe_chrom <- na.omit(seg_chrom_all[,2:16])

#jaccard similarity index
jaccard <- dist(cluster_dataframe_2, method = "binary", upper = TRUE, diag = TRUE)
jaccard_df <- as.data.frame(1-as.matrix(jaccard))
plot(hclust(jaccard))
heatmap(as.matrix(jaccard))

#cluster and dendrogram
autoplot(prcomp(cluster_dataframe_2), label= TRUE, shape= FALSE)
dd <- dist(scale(cluster_dataframe_2), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
plot(hc)


    ##selecting significant samples based on a threshold of (-0.33, +0.33) for CNV 
    ## This selects the samples that have 
#only in either group
cancer_types_sig_threshold <- cancer_types_sig[abs(cancer_types_sig$Mean) > 0.33,]
#in both groups
cancer_types_sig_threshold_both_groups <- data.frame()
for (i in unique(cancer_types_sig$cancer_type.x)){
  cancer_type_quiescent <- cancer_types_sig[cancer_types_sig$cancer_type.x == i & cancer_types_sig$group == "Highly Quiescent",]
  cancer_type_fast.cycling <- cancer_types_sig[cancer_types_sig$cancer_type.x == i & cancer_types_sig$group == "Fast Cycling",]
  for (j in unique(intersect(cancer_type_quiescent$chrom.x, cancer_type_fast.cycling$chrom.x))){
    for (h in unique(intersect(cancer_type_quiescent[cancer_type_quiescent$chrom.x == j,]$segment, cancer_type_fast.cycling[cancer_type_fast.cycling$chrom.x == j,]$segment))){
      if (abs(cancer_type_quiescent[cancer_type_quiescent$chrom.x == j & cancer_type_quiescent$segment == h,]$Mean) >= 0.33 & 
          abs(cancer_type_fast.cycling[cancer_type_fast.cycling$chrom.x == j & cancer_type_fast.cycling$segment == h,]$Mean) >= 0.33){
        cancer_types_sig_threshold_both_groups <- rbind(cancer_types_sig_threshold_both_groups, cancer_type_quiescent[cancer_type_quiescent$chrom.x == j & cancer_type_quiescent$segment == h,], 
                                            cancer_type_fast.cycling[cancer_type_fast.cycling$chrom.x == j & cancer_type_fast.cycling$segment == h,])
      }
    }
  }
}

#matching cancer_type_sig_threshold segments with genes within these segments
gr_cancer_types_sig_threshold <- GRanges(
  seqnames = Rle(cancer_types_sig_threshold$chrom.x),
  ranges = IRanges(cancer_types_sig_threshold$start.pos, end= cancer_types_sig_threshold$end.pos),
  CancerType = cancer_types_sig_threshold$cancer_type.x,
  CNV = cancer_types_sig_threshold$Mean,
  Group = cancer_types_sig_threshold$group)

gr_hg38_genome <- GRanges(
  seqnames = Rle(hg38_genome$Chromosome.scaffold.name),
  ranges = IRanges(hg38_genome$Gene.start..bp., end= hg38_genome$Gene.end..bp.),
  genes = hg38_genome$Gene.name 
)


overlaps <- as.data.frame(findOverlaps(gr_cancer_types_sig_threshold, gr_hg38_genome))
overlaps_genes <- cbind(as.data.frame(gr_hg38_genome[overlaps$subjectHits], row.names = seq(from = 1, to = nrow(overlaps), 1)), as.data.frame(gr_cancer_types_sig_threshold[overlaps$queryHits], row.names = seq(from = 1, to = nrow(overlaps), by = 1)))
overlaps_genes <- overlaps_genes[,6:length(overlaps_genes)]
overlaps_genes <- overlaps_genes[!overlaps_genes$genes == "",]
overlaps_genes$segment <- paste(overlaps_genes$start, "-", overlaps_genes$end)
colnames(overlaps_genes)[2] <- "chrom"

overlaps_genes_chr <- overlaps_genes %>% group_split(overlaps_genes$chrom)
aggregate(data= overlaps_genes_chr[[1]], genes ~ segment, FUN= length)

    ##summary statistics
#finding how many counts each chromosome has
genes_per_chromosome <- data.frame()
for (i in unique(overlaps_genes$chrom)){
  genes_per_chromosome <- rbind(genes_per_chromosome, c(i, length(unique(overlaps_genes[overlaps_genes$chrom == i,]$genes)), 
                                                        length(unique(overlaps_genes[overlaps_genes$chrom == i & overlaps_genes$CNV < 0,]$genes)),     
                                                        length(unique(overlaps_genes[overlaps_genes$chrom == i & overlaps_genes$CNV > 0,]$genes))))
}
colnames(genes_per_chromosome) <- c("Chromosome", "Unique genes", "Deleted Genes", "Amplified Genes")
genes_per_chromosome$Chromosome <- as.numeric(genes_per_chromosome$Chromosome)
#dataframe for plot
genes_per_chromosome_plot_df <- data.frame()
for (i in unique(overlaps_genes$chrom)){
  genes_per_chromosome_plot_df <- rbind(genes_per_chromosome_plot_df, c(i, length(unique(overlaps_genes[overlaps_genes$chrom == i & overlaps_genes$CNV < 0,]$genes)), "Deleted"))
  genes_per_chromosome_plot_df <- rbind(genes_per_chromosome_plot_df, c(i, length(unique(overlaps_genes[overlaps_genes$chrom == i & overlaps_genes$CNV > 0,]$genes)), "Amplified"))
}
colnames(genes_per_chromosome_plot_df) <- c("Chromosome", "Genes", "CNV")

ggplot(genes_per_chromosome_plot_df, aes(x = Chromosome, y = Genes, fill = CNV)) +
        geom_bar(stat= "identity")+ 
        ggtitle("Count of affected genes per chromosome arm")

#finding out which segments are shared between cancers
cancer_types_sig_threshold$segment_unique <- paste(cancer_types_sig_threshold$segment,"_",cancer_types_sig_threshold$chrom.x)
unique_segment_dataframe <- as.data.frame(table(cancer_types_sig_threshold$segment_unique))
table(unique_segment_dataframe$Freq)
#adding the cancers where this segment is found
for (i in 1:nrow(unique_segment_dataframe)){
  unique_segment_dataframe$Cancers[i] <- paste(as.character(as.vector(unique(cancer_types_sig_threshold[cancer_types_sig_threshold$segment_unique == unique_segment_dataframe[i,]$Var1,]$cancer_type.x))), sep="' '", collapse=", ")
  unique_segment_dataframe$CNV[i] <- paste(as.character(as.vector(unique(cancer_types_sig_threshold[cancer_types_sig_threshold$segment_unique == unique_segment_dataframe[i,]$Var1,]$Mean))), sep="' '", collapse=", ")
}



#enrichR analysis
#setting databases
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- c("KEGG_2021_Human", "Reactome_2016", "BioCarta_2016")

#searching for enriched genes overall, in quiescent samples, and in fast cycling samples
if (websiteLive) {
  enriched_overall <- enrichr(overlaps_genes$genes, dbs)
}

if (websiteLive) {
  enriched_quiescent <- enrichr(overlaps_genes[overlaps_genes$Group == "Highly Quiescent",]$genes, dbs)
}

if (websiteLive) {
  enriched_fast.cycling <- enrichr(overlaps_genes[overlaps_genes$Group == "Fast Cycling",]$genes, dbs)
}

#visualising enrichment analysis results and saving plots
pdf(file= "Enrichment_KEGG_Reactome_BioCarta.pdf" )
par( mfrow= c(3,1) )

if (websiteLive) plotEnrich(enriched_overall[["KEGG_2021_Human"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
if (websiteLive) plotEnrich(enriched_quiescent[["KEGG_2021_Human"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
if (websiteLive) plotEnrich(enriched_fast.cycling[["KEGG_2021_Human"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

if (websiteLive) plotEnrich(enriched_overall[["Reactome_2016"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
if (websiteLive) plotEnrich(enriched_quiescent[["Reactome_2016"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
if (websiteLive) plotEnrich(enriched_fast.cycling[["Reactome_2016"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

if (websiteLive) plotEnrich(enriched_overall[["BioCarta_2016"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
if (websiteLive) plotEnrich(enriched_quiescent[["BioCarta_2016"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
if (websiteLive) plotEnrich(enriched_fast.cycling[["BioCarta_2016"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

dev.off()


#separating genes into amplified and deleted genes
#carrying out enrichment analysis
overlaps_genes_deleted <- overlaps_genes[overlaps_genes$CNV < 0,]
overlaps_genes_amplified <- overlaps_genes[overlaps_genes$CNV > 0,]

setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- c("KEGG_2021_Human", "Reactome_2016", "BioCarta_2016")

if (websiteLive) {
  enriched_deleted_quiescent <- enrichr(overlaps_genes_deleted[overlaps_genes_deleted$Group == "Highly Quiescent",]$genes, dbs)
}
if (websiteLive) {
  enriched_deleted_fast.cycling <- enrichr(overlaps_genes_deleted[overlaps_genes_deleted$Group == "Fast Cycling",]$genes, dbs)
}
if (websiteLive) {
  enriched_amplified_quiescent <- enrichr(overlaps_genes_amplified[overlaps_genes_amplified$Group == "Highly Quiescent",]$genes, dbs)
}
if (websiteLive) {
  enriched_amplified_fast.cycling <- enrichr(overlaps_genes_amplified[overlaps_genes_amplified$Group == "Fast Cycling",]$genes, dbs)
}

#plotting KEGG
pdf(file= "Enrichment_KEGG_del_HQ_del_FC_ampl_HQ_ampl_FC.pdf")
if (websiteLive) plotEnrich(enriched_deleted_quiescent[["KEGG_2021_Human"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
if (websiteLive) plotEnrich(enriched_deleted_fast.cycling[["KEGG_2021_Human"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
if (websiteLive) plotEnrich(enriched_amplified_quiescent[["KEGG_2021_Human"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
if (websiteLive) plotEnrich(enriched_amplified_fast.cycling[["KEGG_2021_Human"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

#plotting Reactome
#plotting KEGG
pdf(file= "Enrichment_Reactome_del_HQ_del_FC_ampl_HQ_ampl_FC.pdf")
if (websiteLive) plotEnrich(enriched_deleted_quiescent[["Reactome_2016"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
if (websiteLive) plotEnrich(enriched_deleted_fast.cycling[["Reactome_2016"]], showTerms = 10, numChar = 30, y = "Count", orderBy = "P.value")
if (websiteLive) plotEnrich(enriched_amplified_quiescent[["Reactome_2016"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
if (websiteLive) plotEnrich(enriched_amplified_fast.cycling[["Reactome_2016"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

#plotting BioCarta
#plotting KEGG
pdf(file= "Enrichment_BioCarta_del_HQ_del_FC_ampl_HQ_ampl_FC.pdf")
if (websiteLive) plotEnrich(enriched_deleted_quiescent[["BioCarta_2016"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
if (websiteLive) plotEnrich(enriched_deleted_fast.cycling[["BioCarta_2016"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
if (websiteLive) plotEnrich(enriched_amplified_quiescent[["BioCarta_2016"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
if (websiteLive) plotEnrich(enriched_amplified_fast.cycling[["BioCarta_2016"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

#saving the results in dataframes
enriched_amplified_fast.cycling_KEGG_df <- enriched_amplified_fast.cycling[[1]]
enriched_amplified_quiescent_KEGG_df <- enriched_amplified_quiescent[[1]]
enriched_deleted_fast.cycling_KEGG_df <- enriched_deleted_fast.cycling[[1]]
enriched_deleted_quiescent_KEGG_df <- enriched_deleted_quiescent[[1]]

enriched_amplified_fast.cycling_Reactome_df <- enriched_amplified_fast.cycling[[2]]
enriched_amplified_quiescent_Reactome_df <- enriched_amplified_quiescent[[2]]
enriched_deleted_fast.cycling_Reactome_df <- enriched_deleted_fast.cycling[[2]]
enriched_deleted_quiescent_Reactome_df <- enriched_deleted_quiescent[[2]]

enriched_amplified_fast.cycling_BioCarta_df <- enriched_amplified_fast.cycling[[3]]
enriched_amplified_quiescent_BioCarta_df <- enriched_amplified_quiescent[[3]]
enriched_deleted_fast.cycling_BioCarta_df <- enriched_deleted_fast.cycling[[3]]
enriched_deleted_quiescent_BioCarta_df <- enriched_deleted_quiescent[[3]]

for (i in unique(pancancer_quiescence_groups$CancerType)){
  print(i)
  print(table(pancancer_quiescence_groups[pancancer_quiescence_groups$CancerType == i,]$Group))
}

#checking whether the deletion/ amplification of genes matches with the expression of mRNA (RSEM normalised by Illumina, downloaded from cBioPortal)
setwd("~/Desktop/Fourth year project/Data/CNV/mRNA expression normalised by Illumina/All")
files_mRNA_expression <- list.files(pattern = "txt")

dataframe_mRNA_expression_pancancer <- data.frame()
for (i in 1:length(files_mRNA_expression)){
  file_expression <- read.delim(files_mRNA_expression[i])
  file_expression$SAMPLE_ID <- gsub(" ", "", paste(file_expression$SAMPLE_ID, "A"))
  file_expression_merged_pancancer <- merge(file_expression, pancancer_quiescence_groups, by.x= "SAMPLE_ID", by.y= "Barcode2")
  
  file_expression_long <- melt(file_expression_merged_pancancer[, -which(names(file_expression_merged_pancancer) %in% c("X.1", "X", "Barcode", "STUDY_ID"))], id.vars = c("SAMPLE_ID", "Group", "QuiescenceScore", "CancerType"))
  
  dataframe_mRNA_expression_pancancer <- rbind(file_expression_long, dataframe_mRNA_expression_pancancer)
}
rm(file_expression, file_expression_long, file_expression_merged_pancancer)
dataframe_mRNA_expression_pancancer$value <- log2(dataframe_mRNA_expression_pancancer$value + 1)

#checking how many of the genes in the overlaps genes dataframe were covered from the data in cbioportal
dataframe_genes_queried <- data.frame()
for (i in unique(overlaps_genes$CancerType)){
  dataframe_genes_queried <- rbind(dataframe_genes_queried, 
                                   cbind(CancerType = i,
                                         CNV_N_unique_genes = length(unique(overlaps_genes[overlaps_genes$CancerType == i,]$genes)),
                                         mRNA_N_unique_genes = length(unique(dataframe_mRNA_expression_pancancer[dataframe_mRNA_expression_pancancer$CancerType == i,]$variable))))
}

#removing the NAs from mRNA expression values
dataframe_mRNA_expression_pancancer <- na.omit(dataframe_mRNA_expression_pancancer)

dataframe_genes_queried_after_na_omit <- c()
for (i in unique(overlaps_genes$CancerType)){
  dataframe_genes_queried_after_na_omit <- c(dataframe_genes_queried_after_na_omit, 
                                             length(unique(dataframe_mRNA_expression_pancancer[dataframe_mRNA_expression_pancancer$CancerType == i,]$variable)))
}
dataframe_genes_queried$mRNA_N_unique_genes_NA_omit <- dataframe_genes_queried_after_na_omit

#creating a dataframe with the aggregated values for each variable, by group and by cancer
dataframe_genes_mRNA_expression_aggregated <- dataframe_mRNA_expression_pancancer %>% group_by(variable, CancerType, Group) %>% summarise(mean = mean(value))

#finding which genes have groups of 0 in both groups
genes_with_mean_0 <- data.frame()
for (j in unique(dataframe_genes_mRNA_expression_aggregated$CancerType)){
  for (i in unique(dataframe_genes_mRNA_expression_aggregated[dataframe_genes_mRNA_expression_aggregated$CancerType == j,]$variable)){
    if(dataframe_genes_mRNA_expression_aggregated[dataframe_genes_mRNA_expression_aggregated$variable == i &
                                                  dataframe_genes_mRNA_expression_aggregated$Group == "Highly Quiescent" & 
                                                  dataframe_genes_mRNA_expression_aggregated$CancerType == j,]$mean == 0 & 
       dataframe_genes_mRNA_expression_aggregated[dataframe_genes_mRNA_expression_aggregated$variable == i &
                                                  dataframe_genes_mRNA_expression_aggregated$Group == "Fast Cycling" & 
                                                  dataframe_genes_mRNA_expression_aggregated$CancerType == j,]$mean == 0){
      genes_with_mean_0 <- rbind(genes_with_mean_0, cbind(cancer = j, gene = i))
    }
  }
}
print(genes_with_mean_0)

#removing those genes with mean 0 in both groups for wilcoxon test (t test works fine)
dataframe_mRNA_expression_pancancer <- dataframe_mRNA_expression_pancancer[!dataframe_mRNA_expression_pancancer$variable %in% genes_with_mean_0$gene,]


#testing the expression between groups for each gene (variable)
mRNA_expression_test_dataframe <- dataframe_mRNA_expression_pancancer %>% group_by(variable, CancerType) %>% 
  wilcox_test(value ~ Group, ref.group= "Highly Quiescent") %>% adjust_pvalue(method = "fdr") %>% add_significance(p.col = "p.adj")

#checking whether the direction of the change between CNV and RSEM mRNA expression is the same:
#deletions -> lower expression and amplifications -> higher expression
overlaps_genes_merged_mRNA_expression <- na.omit(merge(overlaps_genes, mRNA_expression_test_dataframe, by.x = c("CancerType", "genes"), by.y = c("CancerType", "variable")))
overlaps_genes_merged_mRNA_expression <- overlaps_genes_merged_mRNA_expression[,-c(6,7,12:15,17)]

for (i in 1:nrow(overlaps_genes_merged_mRNA_expression)){
  if (overlaps_genes_merged_mRNA_expression$p.adj[i] < 0.05){
    if (overlaps_genes_merged_mRNA_expression$statistic[i] > 0 & overlaps_genes_merged_mRNA_expression$CNV[i] > 0){
      overlaps_genes_merged_mRNA_expression$direction[i] <- "same"
    }
    if (overlaps_genes_merged_mRNA_expression$statistic[i] < 0 & overlaps_genes_merged_mRNA_expression$CNV[i] < 0){
      overlaps_genes_merged_mRNA_expression$direction[i] <- "same"
    } 
    else {
      overlaps_genes_merged_mRNA_expression$direction[i] <- "opposite"
    }
  } else {
    overlaps_genes_merged_mRNA_expression$direction[i] <- NA
  }
}

#merging with deleted genes quiescent
quiescent_deleted_overlaps_and_mRNA <-  merge(overlaps_genes[overlaps_genes$Group == "Highly Quiescent" & 
                 overlaps_genes$CNV < 0,], mRNA_expression_test_dataframe, by.x = c("genes", "CancerType"), by.y = c("variable", "CancerType"))

quiescent_deleted_overlaps_and_mRNA$genes_cancer <- paste(quiescent_deleted_overlaps_and_mRNA$genes, "_", quiescent_deleted_overlaps_and_mRNA$CancerType)

#merging with amplified genes quiescent
quiescent_amplified_overlaps_and_mRNA <-  merge(overlaps_genes[overlaps_genes$Group == "Highly Quiescent" & 
                                                               overlaps_genes$CNV > 0,], mRNA_expression_test_dataframe, by.x = c("genes", "CancerType"), by.y = c("variable", "CancerType"))

quiescent_amplified_overlaps_and_mRNA$genes_cancer <- paste(quiescent_amplified_overlaps_and_mRNA$genes, "_", quiescent_amplified_overlaps_and_mRNA$CancerType)


dataframe_mRNA_expression_pancancer$genes_cancer <- paste(dataframe_mRNA_expression_pancancer$variable, "_", dataframe_mRNA_expression_pancancer$CancerType)

#plotting the quiescent deleted genes
mRNA_expression_test_dataframe$genes_cancer <- paste(mRNA_expression_test_dataframe$variable, "_", mRNA_expression_test_dataframe$CancerType)



for (i in unique(quiescent_deleted_overlaps_and_mRNA$genes_cancer)[[3]]){
  print(ggplot(data = dataframe_mRNA_expression_pancancer[dataframe_mRNA_expression_pancancer$genes_cancer == i,], aes(x = Group, y = value)) + 
          geom_boxplot(aes(fill = Group)) +
          ggtitle(quiescent_deleted_overlaps_and_mRNA[quiescent_deleted_overlaps_and_mRNA$genes_cancer == i,]$genes_cancer[[1]]) + geom_jitter(size = 1, aes(color = Group)) +
          stat_pvalue_manual(data = mRNA_expression_test_dataframe[mRNA_expression_test_dataframe$genes_cancer ==i,], label = "p.adj.signif", y.position= 11, size = 9) +
          theme(axis.title=element_text(size=14,face="bold"),
                axis.text=element_text(size=12),
                plot.title = element_text(size = 15, face = "bold"),
                legend.text = element_text(size=12),
                legend.position = "top"))
}

theme <- theme(axis.title=element_text(size=14,face="bold"),
               axis.text=element_text(size=12),
               plot.title = element_text(size = 15, face = "bold"),
               legend.text = element_text(size=12),
               legend.position = "top")

i <- "CDH1 _ BRCA"
p1 <- ggplot(data = dataframe_mRNA_expression_pancancer[dataframe_mRNA_expression_pancancer$genes_cancer == i,], aes(x = Group, y = value)) + 
  geom_boxplot(aes(fill = Group)) +
  ggtitle(quiescent_deleted_overlaps_and_mRNA[quiescent_deleted_overlaps_and_mRNA$genes_cancer == i,]$genes_cancer[[1]]) + geom_jitter(size = 1, aes(color = Group)) +
  stat_pvalue_manual(data = mRNA_expression_test_dataframe[mRNA_expression_test_dataframe$genes_cancer ==i,], label = "p.adj.signif", y.position= 17, size = 8) + theme + 
  ylab("log(z-score)")

i <- "CDH3 _ BRCA"
p2 <- ggplot(data = dataframe_mRNA_expression_pancancer[dataframe_mRNA_expression_pancancer$genes_cancer == i,], aes(x = Group, y = value)) + 
  geom_boxplot(aes(fill = Group)) +
  ggtitle(quiescent_deleted_overlaps_and_mRNA[quiescent_deleted_overlaps_and_mRNA$genes_cancer == i,]$genes_cancer[[1]]) + geom_jitter(size = 1, aes(color = Group)) +
  stat_pvalue_manual(data = mRNA_expression_test_dataframe[mRNA_expression_test_dataframe$genes_cancer ==i,], label = "p.adj.signif", y.position= 15, size = 8) + theme + 
  ylab("log(z-score)")

ggarrange(p1, p2, ncol = 2)









#plotting CNV, SD/segment, SD distribution, and Mean distribution for each chromosome, same cancer
pdf(file= "/Users/victorkonstantellos/Desktop/Fourth year project/Data/CNV/BRCA_CNV_100kbp_chromosome_segments.pdf")
par(mfrow=c(length(unique(overlap_dataframe_aggregated_all$chrom)),4))
for (i in unique(overlap_dataframe_aggregated_all$chrom)){
  print(ggplot(data= overlap_dataframe_aggregated_all[overlap_dataframe_aggregated_all$chrom == i,], mapping = aes(x= start.pos, y= Mean, fill = group))+
          geom_area(position = "identity", alpha = 0.5)+
          scale_fill_brewer(palette="Dark2")+
          theme(legend.position="top")+
          ggtitle(paste("CNV, Chr", i, ", 100kbp segments, BRCA"))+ xlab("Position, 100kbp segments")+ ylab("Mean CNV"))
  
  print(ggplot(data= overlap_dataframe_aggregated_all[overlap_dataframe_aggregated_all$chrom == i,], mapping = aes(x = start.pos, y= SD, fill = group))+
          geom_area(position = "identity", alpha = 0.5)+
          scale_fill_brewer(palette="Dark2")+
          theme(legend.position="top")+
          ggtitle(paste("SD Chr", i, ", 100kbp segments, BRCA"))+ xlab("Position, 100kbp segments")+ ylab("SD"))
  
  print(ggplot(data= overlap_dataframe_aggregated_all[overlap_dataframe_aggregated_all$chrom == i,], mapping = aes(x = SD, fill = group))+
          geom_density(alpha= 0.5)+
          scale_fill_manual(values=c("cyan2", "brown1"))+
          ggtitle(paste("Distribution of SD Chr", i, ", BRCA")))
  
  print(ggplot(data= overlap_dataframe_aggregated_all[overlap_dataframe_aggregated_all$chrom == i,], mapping = aes(x = Mean, fill = group))+
          geom_density(alpha= 0.5)+
          scale_fill_manual(values=c("cyan2", "brown1"))+
          ggtitle(paste("Distribution of Mean, Chr", i, ", BRCA")))
}
dev.off()

#plotting significant segments only
pdf(file= "/Users/victorkonstantellos/Desktop/Fourth year project/Data/CNV/LUAD_100kbp_singificant_segments.pdf")
par(mfrow=c(length(unique(luad_aggregated_significant$chrom)),1))
for (i in unique(luad_aggregated_significant$chrom)){
  print(ggplot(data= luad_aggregated_significant[luad_aggregated_significant$chrom == i,], mapping = aes(x= start.pos, y= Mean, fill = group))+
          geom_area(position = "identity", alpha = 0.5)+
          scale_fill_brewer(palette="Dark2")+
          theme(legend.position="top")+
          ggtitle(paste("CNV, Chr", i, "Significant segments only, LUAD"))+ xlab("Position, 100kbp segments")+ ylab("Mean CNV"))
}
dev.off()

plot <- ggplot(data = chr5, mapping = aes(x = start.pos, y = Mean, group = group)) + ylim(-0.2, 0.2)+
  geom_area(data = subset(chr5, cancer_type == "BRCA"), aes(color = cancer_type, fill = group), alpha = 0.3, position = "identity") + 
  geom_area(data = subset(chr5, cancer_type == "LUAD"), aes(color = cancer_type, fill = group), alpha = 0.5, position = "identity") + 
  scale_fill_brewer(palette="Dark2")+
  scale_colour_manual(values=c("red", "blue"))+
  theme(legend.position="top")+
  ggtitle("Chr 5, Significant positions, BRCA and LUAD")
for (i in 1:nrow(brca_luad_fc_hq_test_sig_chr5)){
  plot <- plot + annotate("segment", x = brca_luad_fc_hq_test_sig_chr5$start.pos[i], xend = brca_luad_fc_hq_test_sig_chr5$end.pos[i], y = 0.175, yend = 0.175, colour = "blue")
}
print(plot)


for (i in unique(overlaps_genes$CancerType)){
  for (j in unique(overlaps_genes[overlaps_genes$CancerType == i,]$chrom)){
    print(c(length(unique(overlaps_genes[overlaps_genes$CancerType == i & overlaps_genes$chrom == j & overlaps_genes$CNV < 0,]$segment)), i, j))
    print(c(length(unique(overlaps_genes[overlaps_genes$CancerType == i & overlaps_genes$chrom == j & overlaps_genes$CNV > 0,]$segment)), i, j))
  }
}

#summary statistics
blank_theme <- theme_minimal()+
  theme(
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=16, face="bold"),
    legend.title=element_text(size=15),
    legend.text=element_text(size=15),
    legend.position = "top",
  )

cbPalette <- brewer.pal(11, "Paired")

cancer_types_sig_threshold$Ampl_or_Del <- ifelse(cancer_types_sig_threshold$Mean > 0, yes = "Amplification", no = "Deletion")

summary <- as.data.frame(rbind(as.data.frame(table(cancer_types_sig_threshold[cancer_types_sig_threshold$Ampl_or_Del == "Amplification",]$group, 
                                                   cancer_types_sig_threshold[cancer_types_sig_threshold$Ampl_or_Del == "Amplification",]$cancer_type.x, 
                                                   cancer_types_sig_threshold[cancer_types_sig_threshold$Ampl_or_Del == "Amplification",]$chrom.x)), 
                               as.data.frame(table(cancer_types_sig_threshold[cancer_types_sig_threshold$Ampl_or_Del == "Deletion",]$group, 
                                                   cancer_types_sig_threshold[cancer_types_sig_threshold$Ampl_or_Del == "Deletion",]$cancer_type.x, 
                                                   cancer_types_sig_threshold[cancer_types_sig_threshold$Ampl_or_Del == "Deletion",]$chrom.x))))
colnames(summary)[2:4] <- c("Cancer", "Chr", "Segments")
summary$Chr <- factor(summary$Chr, levels = c("1", "3", "4", "6", "8", "9", "10", "11", "12", "16", "17"))

summary_1 <- as.data.frame(table(cancer_types_sig_threshold$chrom.x))
colnames(summary_1) <- c("Chr", "Segments")
summary_1$Chr <- factor(summary_1$Chr, levels = c("1", "3", "4", "6", "8", "9", "10", "11", "12", "16", "17"))


summary_2 <- as.data.frame(table(cancer_types_sig_threshold$group, cancer_types_sig_threshold$Ampl_or_Del))
colnames(summary_2)[2:3] <- c("CNV", "Segments")

summary_3 <- as.data.frame(table(cancer_types_sig_threshold[cancer_types_sig_threshold$group == "Fast Cycling",]$Ampl_or_Del, cancer_types_sig_threshold[cancer_types_sig_threshold$group == "Fast Cycling",]$chrom.x))
colnames(summary_3) <- c("CNV", "Chr", "Segments")

summary_4 <- as.data.frame(table(cancer_types_sig_threshold[cancer_types_sig_threshold$group == "Highly Quiescent",]$Ampl_or_Del, cancer_types_sig_threshold[cancer_types_sig_threshold$group == "Highly Quiescent",]$chrom.x))
colnames(summary_4) <- c("CNV", "Chr", "Segments")

p1 <- ggplot(summary_2[summary_2$Var1 == "Fast Cycling",], aes(x="", y=Segments, fill = CNV))+geom_bar(width = 1, stat = "identity") + 
  coord_polar(theta = "y") + blank_theme +  geom_text(aes(label = Segments), position = position_stack(vjust = 0.5), size = 6) + 
  ggtitle("A) Fast Cycling Tumours") 

p2 <- ggplot(summary_2[summary_2$Var1 == "Highly Quiescent",], aes(x="", y=Segments, fill = CNV))+geom_bar(width = 1, stat = "identity") + 
  coord_polar(theta = "y") + blank_theme +  geom_text(aes(label = Segments), position = position_stack(vjust = 0.5), size = 6) + 
  ggtitle("B) Highly Quiescent Tumours") 

p3 <- ggplot(summary[summary$Var1 == "Fast Cycling",], aes(fill= Chr, y= Segments, x= Cancer)) + 
  geom_bar(position="stack", stat="identity") + ggtitle("D) Segment CNV across chromosomes and cancer, Fast Cycling") + 
  theme(plot.title=element_text(size=16, face="bold"),
legend.title=element_text(size=15),
legend.text=element_text(size=15), 
axis.text = element_text(size = 15)) +
  scale_fill_manual(values=cbPalette)

p4 <- ggplot(summary[summary$Var1 == "Highly Quiescent",], aes(fill= Chr, y= Segments, x= Cancer)) + 
  geom_bar(position="stack", stat="identity") + ggtitle("E) Segment CNV across chromosomes and cancer, Highly Quiescent") +
  theme(plot.title=element_text(size=16, face="bold"),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15), 
        axis.text = element_text(size = 15)) + 
  scale_fill_manual(values=cbPalette)

p5 <- ggplot(summary_1, aes(x="", y=Segments, fill = Chr))+geom_bar(width = 1, stat = "identity") + 
  coord_polar(theta = "y") + blank_theme  + theme(legend.position = "right") +
  ggtitle("C) Distibution of CNV in chromosomes") + 
  scale_fill_manual(values=cbPalette) 



grid.arrange(arrangeGrob(p1, p2, p5, ncol = 3),  arrangeGrob(p3, p4, ncol = 2), nrow = 2)


#enrichment plots
m1 <- plotEnrich(enriched_deleted_fast.cycling[["Reactome_2016"]], showTerms = 10, numChar = 30, y = "Count", 
                 orderBy = "P.value", title = "Deleted Genes, FC, Reactome")
m1[["data"]]$P.value <- m1[["data"]]$Adjusted.P.value
m1[["theme"]][["axis.text"]]$size <- 16
m1[["theme"]][["plot.title"]]$size <- 16
m1[["theme"]][["plot.title"]]$face <- "bold"
m1[["theme"]][["axis.text"]]$size <- 15
m1[["theme"]][["legend.text"]]$size <- 15
m1[["theme"]][["legend.title"]]$size <- 15

m2 <- plotEnrich(enriched_deleted_quiescent[["Reactome_2016"]], showTerms = 10, numChar = 30, y = "Count", 
                 orderBy = "P.value", title = "Deleted Genes, HQ, Reactome")
m2[["data"]]$P.value <- m2[["data"]]$Adjusted.P.value
m2[["theme"]][["axis.text"]]$size <- 16
m2[["theme"]][["plot.title"]]$size <- 16
m2[["theme"]][["plot.title"]]$face <- "bold"
m2[["theme"]][["axis.text"]]$size <- 15
m2[["theme"]][["legend.text"]]$size <- 15
m2[["theme"]][["legend.title"]]$size <- 15


m3 <- plotEnrich(enriched_deleted_fast.cycling[["KEGG_2021_Human"]], showTerms = 10, numChar = 30, y = "Count", 
                 orderBy = "P.value", title = "Deleted Genes, FC, KEGG")
m3[["data"]]$P.value <- m3[["data"]]$Adjusted.P.value
m3[["theme"]][["axis.text"]]$size <- 16
m3[["theme"]][["plot.title"]]$size <- 16
m3[["theme"]][["plot.title"]]$face <- "bold"
m3[["theme"]][["axis.text"]]$size <- 15
m3[["theme"]][["legend.text"]]$size <- 15
m3[["theme"]][["legend.title"]]$size <- 15



