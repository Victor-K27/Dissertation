#reading the files and loading libraries
library(readxl)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggplot2)
library(gridExtra)

#reading the files
setwd("~/Desktop/Fourth year project/Data/BRCA1, BRCA2")
brca <- read_excel("BRCA.xlsx", sheet= 1)
pancancer_quiescence_groups <- read.csv2("PancancerQuiescenceGroupsCorrected.csv")

#merging the brca and pancancer dataframes
brca_merged_pancancer <- merge(brca, pancancer_quiescence_groups, by.x = "sample", by.y = "Barcode2")

#making separate dataframes for each condition in the BRCA1 germline 
brca1_somatic_null <- data.frame(brca_merged_pancancer$sample, brca_merged_pancancer$BRCA1_somatic_null, brca_merged_pancancer$QuiescenceScore)
colnames(brca1_somatic_null) <- c("sample", "Condition", "Quiescence_score")
brca1_somatic_null <- brca1_somatic_null[!brca1_somatic_null$Condition == "NA",]
brca1_somatic_null <- brca1_somatic_null[!brca1_somatic_null$Condition == "0",]
brca1_somatic_null <- brca1_somatic_null[!brca1_somatic_null$Condition == "NaN",]

brca1_germ_bi_allelic <- data.frame(brca_merged_pancancer$sample, brca_merged_pancancer$BRCA1_germ_bi_allelic, brca_merged_pancancer$QuiescenceScore)
colnames(brca1_germ_bi_allelic) <- c("sample", "Condition", "Quiescence_score")
brca1_germ_bi_allelic <- brca1_germ_bi_allelic[!brca1_germ_bi_allelic$Condition == "NA",]
brca1_germ_bi_allelic <- brca1_germ_bi_allelic[!brca1_germ_bi_allelic$Condition == "0",]
brca1_germ_bi_allelic <- brca1_germ_bi_allelic[!brca1_germ_bi_allelic$Condition == "NaN",]

brca1_germ_mono_allelic <- data.frame(brca_merged_pancancer$sample, brca_merged_pancancer$BRCA1_germ_mono_allelic, brca_merged_pancancer$QuiescenceScore)
colnames(brca1_germ_mono_allelic) <- c("sample", "Condition", "Quiescence_score")
brca1_germ_mono_allelic <- brca1_germ_mono_allelic[!brca1_germ_mono_allelic$Condition == "NA",]
brca1_germ_mono_allelic <- brca1_germ_mono_allelic[!brca1_germ_mono_allelic$Condition == "0",]
brca1_germ_mono_allelic <- brca1_germ_mono_allelic[!brca1_germ_mono_allelic$Condition == "NaN",]

brca1_deletion <- data.frame(brca_merged_pancancer$sample, brca_merged_pancancer$BRCA1_deletion, brca_merged_pancancer$QuiescenceScore)
colnames(brca1_deletion) <- c("sample", "Condition", "Quiescence_score")
brca1_deletion <- brca1_deletion[!brca1_deletion$Condition == "NA",]
brca1_deletion <- brca1_deletion[!brca1_deletion$Condition == "0",]
brca1_deletion <- brca1_deletion[!brca1_deletion$Condition == "NaN",]

brca1_epigenetic_silencing <- data.frame(brca_merged_pancancer$sample, brca_merged_pancancer$BRCA1_epigenetic_silencing, brca_merged_pancancer$QuiescenceScore)
colnames(brca1_epigenetic_silencing) <- c("sample", "Condition", "Quiescence_score")
brca1_epigenetic_silencing <- brca1_epigenetic_silencing[!brca1_epigenetic_silencing$Condition == "NA",]
brca1_epigenetic_silencing <- brca1_epigenetic_silencing[!brca1_epigenetic_silencing$Condition == "0",]
brca1_epigenetic_silencing <- brca1_epigenetic_silencing[!brca1_epigenetic_silencing$Condition == "NaN",]

brca1_mRNA <- data.frame(brca_merged_pancancer$sample, brca_merged_pancancer$BRCA1_mRNA, brca_merged_pancancer$QuiescenceScore)
colnames(brca1_mRNA) <- c("sample", "Condition", "Quiescence_score")
brca1_mRNA <- brca1_mRNA[!brca1_mRNA$Condition == "NA",]
brca1_mRNA <- brca1_mRNA[!brca1_mRNA$Condition == "0",]
brca1_mRNA <- brca1_mRNA[!brca1_mRNA$Condition == "NaN",]

#grouping somatic
brca1_somatic_grouped <- brca1_somatic_null
brca1_somatic_grouped$Condition <- rep("somatic", nrow(brca1_somatic_grouped))
#grouping epigenetic silencing and mRNA downregulation together
brca1_epigenetic.silencing_mRNA.downregulation <- rbind(brca1_epigenetic_silencing, brca1_mRNA)
brca1_epigenetic.silencing_mRNA.downregulation$Condition <- rep("expression_downregulation", nrow(brca1_epigenetic.silencing_mRNA.downregulation))

#making the brca1 all conditions dataframe
brca1_all.conditions <- rbind (brca1_somatic_grouped, brca1_germ_bi_allelic, brca1_germ_mono_allelic, brca1_deletion, brca1_epigenetic.silencing_mRNA.downregulation)

#adding a column with somatic or germline information
brca1_all.conditions$Somatic_or_germline <- rep(NA, nrow(brca1_all.conditions))

for (i in 1:nrow(brca1_all.conditions)){ifelse (brca1_all.conditions$Condition[i] == "somatic", 
                                                yes= brca1_all.conditions$Somatic_or_germline[i] <- "somatic", 
                                                no = brca1_all.conditions$Somatic_or_germline[i] <- "germline")}

#making separate dataframes for each condition in the BRCA2 germline 
brca2_somatic_null <- data.frame(brca_merged_pancancer$sample, brca_merged_pancancer$BRCA2_somatic_null, brca_merged_pancancer$QuiescenceScore)
colnames(brca2_somatic_null) <- c("sample", "Condition", "Quiescence_score")
brca2_somatic_null <- brca2_somatic_null[!brca2_somatic_null$Condition == "NA",]
brca2_somatic_null <- brca2_somatic_null[!brca2_somatic_null$Condition == "0",]
brca2_somatic_null <- brca2_somatic_null[!brca2_somatic_null$Condition == "NaN",]

brca2_germ_bi_allelic <- data.frame(brca_merged_pancancer$sample, brca_merged_pancancer$BRCA2_germ_bi_allelic, brca_merged_pancancer$QuiescenceScore)
colnames(brca2_germ_bi_allelic) <- c("sample", "Condition", "Quiescence_score")
brca2_germ_bi_allelic <- brca2_germ_bi_allelic[!brca2_germ_bi_allelic$Condition == "NA",]
brca2_germ_bi_allelic <- brca2_germ_bi_allelic[!brca2_germ_bi_allelic$Condition == "0",]
brca2_germ_bi_allelic <- brca2_germ_bi_allelic[!brca2_germ_bi_allelic$Condition == "NaN",]

brca2_germ_undetermined <- data.frame(brca_merged_pancancer$sample, brca_merged_pancancer$BRCA2_germ_undetermined, brca_merged_pancancer$QuiescenceScore)
colnames(brca2_germ_undetermined) <- c("sample", "Condition", "Quiescence_score")
brca2_germ_undetermined <- brca2_germ_undetermined[!brca2_germ_undetermined$Condition == "NA",]
brca2_germ_undetermined <- brca2_germ_undetermined[!brca2_germ_undetermined$Condition == "0",]
brca2_germ_undetermined <- brca2_germ_undetermined[!brca2_germ_undetermined$Condition == "NaN",]

brca2_germ_mono_allelic <- data.frame(brca_merged_pancancer$sample, brca_merged_pancancer$BRCA2_germ_mono_allelic, brca_merged_pancancer$QuiescenceScore)
colnames(brca2_germ_mono_allelic) <- c("sample", "Condition", "Quiescence_score")
brca2_germ_mono_allelic <- brca2_germ_mono_allelic[!brca2_germ_mono_allelic$Condition == "NA",]
brca2_germ_mono_allelic <- brca2_germ_mono_allelic[!brca2_germ_mono_allelic$Condition == "0",]
brca2_germ_mono_allelic <- brca2_germ_mono_allelic[!brca2_germ_mono_allelic$Condition == "NaN",]

brca2_deletion <- data.frame(brca_merged_pancancer$sample, brca_merged_pancancer$BRCA2_deletion, brca_merged_pancancer$QuiescenceScore)
colnames(brca2_deletion) <- c("sample", "Condition", "Quiescence_score")
brca2_deletion <- brca2_deletion[!brca2_deletion$Condition == "NA",]
brca2_deletion <- brca2_deletion[!brca2_deletion$Condition == "0",]
brca2_deletion <- brca2_deletion[!brca2_deletion$Condition == "NaN",]

#grouping  somatic
brca2_somatic_grouped <- brca2_somatic_null
brca2_somatic_grouped$Condition <- rep("somatic", nrow(brca2_somatic_grouped))

#separating brca2 mono allelic, brca2 bi allelic, and undetermined
brca2_germ_mono_allelic$Condition <- rep("mono_germline.null_or_pathogenic", nrow(brca2_germ_mono_allelic))
brca2_germ_bi_allelic$Condition <- rep("bi_germline.null_or_pathogenic", nrow(brca2_germ_bi_allelic))
brca2_germ_undetermined$Condition <- rep("undetermined_germline.null_or_pathogenic", nrow(brca2_germ_undetermined))

#making the brca2 all conditions dataframe
brca2_all.conditions <- rbind(brca2_somatic_grouped, brca2_germ_bi_allelic, brca2_germ_mono_allelic, brca2_germ_undetermined, brca2_deletion)

#adding a column with somatic or germline information
for (i in 1:nrow(brca2_all.conditions)){ifelse (brca2_all.conditions$Condition[i] == "somatic", 
                                                yes= brca2_all.conditions$Somatic_or_germline[i] <- "somatic", 
                                                no = brca2_all.conditions$Somatic_or_germline[i] <- "germline")}

#adding a column with germline for brca1 and brca2
brca1_all.conditions["Germline"] <- rep("brca1", nrow(brca1_all.conditions))
brca2_all.conditions["Germline"] <- rep("brca2", nrow(brca2_all.conditions))

#combining the brca1 and brca2 all conditions in a single dataframe
brca1_brca2_all.conditions <- rbind(brca1_all.conditions, brca2_all.conditions)

#selecting the unmutated from the brca_merged_pancancer
brca_unmutated <- data.frame(sample = brca_merged_pancancer$sample, brca_merged_pancancer[,7:17], QuiescenceScore = brca_merged_pancancer$QuiescenceScore)
test = c()
for (i in 1:nrow(brca_unmutated)) {
  if (all(brca_unmutated[i, 2:12]== 0)) {
    test = append(test, TRUE)
  } else {
    test = append(test, FALSE)
  }
}
brca_unmutated <- brca_unmutated[test,]
brca_unmutated <- brca_unmutated[,c(1,2,13)]
colnames(brca_unmutated) <- c("sample", "Condition", "Quiescence_score")
brca_unmutated[,2] <- rep("wildtype", nrow(brca_unmutated))
brca_unmutated["Germline"] <- rep("wildtype", nrow(brca_unmutated))
brca_unmutated$Somatic_or_germline <- rep("wildtype", nrow(brca_unmutated))
                                  

#binding the brca1/2_all conditions and brca unmutated
brca_all.data <- rbind(brca1_all.conditions, brca2_all.conditions, brca_unmutated)
for (i in 1:nrow(brca_all.data)){ifelse(brca_all.data$Quiescence_score[i] < 0, 
                                        yes = brca_all.data$Group[i] <- "Fast Cycling", 
                                        no = brca_all.data$Group[i] <- "Highly Quiescent")}

#t.test for Quiescence ~ BRCA1_all conditions
brca1_all.conditions <- rbind(brca1_all.conditions, brca_unmutated)
res.t_test.brca1_all.conditions <- brca1_all.conditions %>% t_test(Quiescence_score  ~ Condition, p.adjust.method = "fdr")
print(res.t_test.brca1_all.conditions)
res.anova.brca1_all.conditions <- aov(Quiescence_score ~ Condition, data = brca1_all.conditions)
summary(res.anova.brca1_all.conditions)
ggplot(data = brca1_all.conditions, aes(x = Condition, y = Quiescence_score)) + geom_boxplot() + geom_jitter() +
  ggtitle("Quiescence ~ BRCA1 All Alterations") + stat_pvalue_manual(data = res.t_test.brca1_all.conditions, hide.ns = TRUE, label = "p.adj.signif", y.position = c(11, 13,15)) +
  xlab("Alteration") + ylab("Quiescence Score")

#t.test for Quiescence ~ BRCA2_all conditions
brca2_all.conditions <- rbind(brca2_all.conditions, brca_unmutated)
res.t_test.brca2_all.conditions <- brca2_all.conditions %>% t_test(Quiescence_score  ~ Condition, p.adjust.method = "fdr")
print(res.t_test.brca2_all.conditions)
res.anova.brca2_all.conditions <- aov(Quiescence_score ~ Condition, data = brca2_all.conditions)
summary(res.anova.brca2_all.conditions)
ggplot(data = brca2_all.conditions, aes(x = Condition,  y = Quiescence_score)) + geom_boxplot() + geom_jitter() +
  ggtitle("Quiescence ~ BRCA2 All Alterations") + stat_pvalue_manual(data = res.t_test.brca2_all.conditions, hide.ns = TRUE, label = "p.adj", y.position = c(7)) +
  xlab("Alteration") + ylab("Quiescence Score")


#t.test for Quiescence ~ Germlines
res.t_test.brca_all.data_germline <- brca_all.data %>% t_test(Quiescence_score ~ Germline, p.adjust.method = "fdr")
print(res.t_test.brca_all.data_germline)
res.anova.brca_all.data_germline <- aov(Quiescence_score ~ Germline, data = brca_all.data)
summary(res.anova.brca_all.data_germline)
p1 <- ggplot(data = brca_all.data, aes(x = Germline, y = Quiescence_score)) + geom_boxplot() + geom_jitter() + ggtitle("Quiescence ~ BRCA1/2 Mutations") +
  stat_pvalue_manual(res.t_test.brca_all.data_germline, label = "p.adj.signif", hide.ns = FALSE, y.position = c(9,11,13)) + xlab("BRCA1 and BRCA2 Mutations") + ylab("Quiescence Score")+
  theme(plot.title = element_text(size = 15, face = "bold"), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))

#t.test for Quiescence ~ Type of alteration in BRCA1/2
res.t_test.brca_all.data_condition <- brca_all.data %>% t_test(Quiescence_score ~ Condition, p.adjust.method = "fdr")
print(res.t_test.brca_all.data_condition)
res.anova.brca_all.data_germline <- aov(Quiescence_score ~ Condition, data = brca_all.data)
summary(res.anova.brca_all.data_germline)
ggplot(data = brca_all.data, aes(x = Condition, y = Quiescence_score)) + geom_boxplot() + geom_jitter(aes(color = Germline)) + ggtitle("Quiescence ~ Type of Alteration in BRCA1/2") +
  stat_pvalue_manual(res.t_test.brca_all.data_condition, label = "p.adj.signif", hide.ns = TRUE, y.position = c(11, 13, 15, 17, 19, 21)) + xlab("Type of Alteration") + ylab("Quiescence Score") +
  theme(plot.title = element_text(size = 15, face = "bold"), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))

#t.test for Quiescence ~ WT/ any other condition
#adding a column to brca_all data with WT/mutated information
brca_all.data["WT.or.mutated"] <- rep(NA, nrow(brca_all.data))
for (i in 1:nrow(brca_all.data)){
  if (brca_all.data$Germline[i] != "wildtype"){
    brca_all.data$WT.or.mutated[i] <- "mutated"
  } else {
    brca_all.data$WT.or.mutated[i] <- "wildtype"
  }
}
res.anova.brca_all.data_wt.or.mutated <- aov(Quiescence_score ~ WT.or.mutated, data = brca_all.data)
summary(res.anova.brca_all.data_wt.or.mutated)
res.t_test.brca_all.data_wt.or.mutated <- brca_all.data %>% t_test(Quiescence_score ~ WT.or.mutated, p.adjust.method = "fdr") %>% add_significance(p.col = "p")
print(res.t_test.brca_all.data_wt.or.mutated)
p2 <- ggplot(data = brca_all.data, aes(x = WT.or.mutated, y = Quiescence_score)) + geom_boxplot() + geom_jitter(aes(color= Germline)) + ggtitle("Quiescence ~ WT or Mutated BRCA1/2") +
  stat_pvalue_manual(res.t_test.brca_all.data_wt.or.mutated, label = "p.signif", hide.ns = TRUE, y.position = c(11)) + xlab("WT or Mutated BRCA1/2") + ylab("Quiescence Score") +
  theme(plot.title = element_text(size = 15, face = "bold"), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), legend.position = "top")

#t.test for Quiescence ~ Somatic_or_germline
res.t_test.brca_all.data_somatic.or.germline <- brca_all.data %>% t_test(Quiescence_score ~ Somatic_or_germline, p.adjust.method = "fdr")
print(res.t_test.brca_all.data_somatic.or.germline)
ggplot(data = brca_all.data, aes(x = Somatic_or_germline, y = Quiescence_score)) + geom_boxplot() + geom_jitter(aes(color = Germline)) + ggtitle("Quiescence ~ Somatic or Germline Mutations in BRCA1/2") +
  stat_pvalue_manual(res.t_test.brca_all.data_somatic.or.germline, label = "p.adj.signif", hide.ns = FALSE, y.position = c(11,13,15)) + xlab("Somatic or Germline Mutations in BRCA1/2") + ylab("Quiescence Score") +
  theme(plot.title = element_text(size = 15, face = "bold"), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), legend.position = "top")

res.t_test.brca_all.data_somatic.or.germline <- brca_all.data %>% t_test(Quiescence_score ~ new, p.adjust.method = "fdr")
ggplot(data = brca_all.data, aes(x = new, y = Quiescence_score)) + geom_boxplot() + geom_jitter(aes(color = Germline)) + ggtitle("Quiescence ~ Somatic or Germline Mutations in BRCA1/2") +
  stat_pvalue_manual(res.t_test.brca_all.data_somatic.or.germline, label = "p.adj.signif", hide.ns = TRUE, y.position = 11) + xlab("Somatic or Germline Mutations in BRCA1/2") + ylab("Quiescence Score") +
  theme(plot.title = element_text(size = 15, face = "bold"), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), legend.position = "top")


grid.arrange(p1, p2, nrow = 2)
grid.arrange(p2, p1, nrow = 1)



                            #### -----  Homologous Repair Deficiency Analysis  ----- ####


#loading libraries
library(readxl)
library(ggpubr)
library(stats)
library(ggplot2)

#reading files
setwd("~/Desktop/Fourth year project/Data/HRD")
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
res.t_test.hrd <- t_test(data = brca_merged_pancancer, QuiescenceScore ~ HRD)
print(res.t_test.hrd)
ggplot(data = brca_merged_pancancer, aes(x = HRD, y= QuiescenceScore)) + geom_boxplot() + geom_jitter() + 
  ggtitle("Quiescence ~ Homologous Repair Deficiency") + xlab("Homologous Repair Deficiencu (HRD)") + ylab("Quiescence Score") +
  stat_pvalue_manual(data = res.t_test.hrd, hide.ns = FALSE, y.position = 10, label = "p")

#counting quiescent samples with and without hrd
print(table(brca_merged_pancancer$Group, brca_merged_pancancer$HRD))
count_dataframe <- data.frame(matrix(c(28, 19, 37, 6), nrow= 2, ncol = 2))
colnames(count_dataframe) <- c("FALSE", "TRUE")
rownames(count_dataframe) <- c("Fast Cycling", "Highly Quiescent")
fisher.test(count_dataframe)

#making the plots 
p1 <- ggplot(data = brca_all.data, aes(x = WT.or.mutated, y = Quiescence_score)) + geom_boxplot() + geom_jitter(aes(color= Germline)) + ggtitle("Quiescence ~ WT or Altered BRCA1/2") +
  stat_pvalue_manual(res.t_test.brca_all.data_wt.or.mutated, label = "p", hide.ns = TRUE, y.position = c(11)) + xlab("WT or Altered BRCA1/2") + ylab("Quiescence Score") +
  theme(legend.position = "top")

p2 <- ggplot(data = brca_all.data, aes(x = Germline, y = Quiescence_score)) + geom_boxplot() + geom_jitter() + ggtitle("Quiescence ~ BRCA1 and BRCA2 Alterations") +
  stat_pvalue_manual(res.t_test.brca_all.data_germline, label = "p.adj.signif", hide.ns = FALSE, y.position = c(9,11,13)) + xlab("BRCA1 and BRCA2 Alterations") + ylab("Quiescence Score")

grid.arrange(p1, p2, nrow = 1)


p3 <- ggplot(data = brca1_all.conditions, aes(x = Condition, y = Quiescence_score)) + geom_boxplot() + geom_jitter() +
  ggtitle("Quiescence ~ BRCA1 All Alterations") + stat_pvalue_manual(data = res.t_test.brca1_all.conditions, hide.ns = TRUE, label = "p.adj.signif", y.position = c(11, 13,15)) +
  xlab("Alteration") + ylab("Quiescence Score")

p4 <- ggplot(data = brca2_all.conditions, aes(x = Condition,  y = Quiescence_score)) + geom_boxplot() + geom_jitter() +
  ggtitle("Quiescence ~ BRCA2 All Alterations") + stat_pvalue_manual(data = res.t_test.brca2_all.conditions, hide.ns = TRUE, label = "p.adj", y.position = c(7)) +
  xlab("Alteration") + ylab("Quiescence Score")

grid.arrange(p3, p4, nrow = 2)


