#loading libraries
library(readxl)
library(ggpubr)
library(stats)
library(ggplot2)

#reading files
setwd("C:/Users/victo/Desktop/Fourth year project/Data/HRD")
stad <- read_excel("STAD.xlsx")
brca <- read_excel("BRCA.xlsx", sheet= 3)
pancancer_quiescence_groups <- read.csv2("PancancerQuiescenceGroupsCorrected.csv")

#merging brca and pancancer
brca_merged_pancancer <- merge(brca, pancancer_quiescence_groups, by.x = "sample", by.y = "Barcode2")

#creating binary columns for HRD and Quiescence
for (i in 1:nrow(brca_merged_pancancer)){
  ifelse(test = brca_merged_pancancer$Group[i] == "Highly Quiescent", 
         yes = brca_merged_pancancer$Binary_group[i] <- 1, 
         no = brca_merged_pancancer$Binary_group[i] <- 0)
  
  ifelse(test = brca_merged_pancancer$HRD[i] == "TRUE",
         yes = brca_merged_pancancer$HRD_binary[i] <- 1,
         no = brca_merged_pancancer$HRD_binary[i] <- 0)
}

#correlation test for Group_binary ~ HRD
cor.test(x= as.numeric(brca_merged_pancancer$HRD_binary), y= as.numeric(brca_merged_pancancer$Binary_group))

#correlation test for Quiescence score ~ HRD
cor.test(x= as.numeric(brca_merged_pancancer$HRD_binary), y= as.numeric(brca_merged_pancancer$QuiescenceScore))

#t.test for Quiescence score ~ HRD
t.test(x= as.logical(brca_merged_pancancer$HRD), y= brca_merged_pancancer$QuiescenceScore)


#boxplot for Quiescence ~ HRD
ggboxplot(data= na.omit(brca_merged_pancancer), x= "HRD", y= "QuiescenceScore", add = "jitter", title = "Quiescence ~ HRD, BRCA", color = "HRD", palette = c("red", "blue")) +
  stat_compare_means(comparisons = na.omit(brca), aes("HRD", "QuiescenceScore")) + stat_compare_means(label.y = 10) +
  stat_compare_means(data= na.omit(brca_merged_pancancer), method = "t.test", label.y = 11)


