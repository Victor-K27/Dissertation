#loading the data from the script
setwd("~/Desktop/Fourth year project/Data/CNV/aneuploidy analysis")
load('myEnvironment.RData')

#loading packages
library(readxl)
library(tidyverse)
library(rstatix)
library(ggplot2)
library(ggfortify)
library(reshape)
library(pheatmap)
library(aod)
library(dplyr)
library(pscl)
library(randomForest)
library(caret)
library(adabag)
library(gridExtra)
library(ggpubr)
library(stringr)

#reading pancancer file
setwd("~/Desktop/Fourth year project/Data/CNV")
pancancer_quiescence_groups <- read.csv2("PancancerQuiescenceGroupsCorrected.csv")

#reading aneuploidy file
setwd("~/Desktop/Fourth year project/Data/CNV/aneuploidy analysis")
aneuploidy <- read_excel("TCGA_aneuploidy.xlsx")
aneuploidy$Sample <- gsub(" ", "", paste(aneuploidy$Sample, "A"))
#selecting only whole chromosome arm data
aneuploidy <- aneuploidy[,c(1,2,14:52)]

#merging with pancancer
aneuploidy_merged_pancancer <- merge(aneuploidy, pancancer_quiescence_groups, by.x = "Sample", by.y = "Barcode2")

#creating a dataframe with the balance of HQ and FC samples in all cancer types
balance_all_cancers <- as.data.frame(table(aneuploidy_merged_pancancer$Type, aneuploidy_merged_pancancer$Group))
colnames(balance_all_cancers) <- c("Cancer", "Group", "Freq")

pdf(file= "balance_aneuploidy_pancancer.pdf")
for (i in unique(balance_all_cancers$Cancer)){
  print(ggplot(balance_all_cancers[balance_all_cancers$Cancer == i,], aes(x= "", y= Freq, fill= Group))+
    geom_col(color= "black")+
    geom_text(aes(label= Freq), position= position_stack(vjust = 0.5))+
    coord_polar(theta= "y")+
    theme_void()+
    ggtitle(i))
}
dev.off()

#removing cancers that have less than 3 samples in either group
aneuploidy_merged_pancancer <- aneuploidy_merged_pancancer[!aneuploidy_merged_pancancer$Type %in% c("GBM", "TGCT", "THYM", "CHOL", "UCS", "ACC", "PCPG"),]

#splitting the aneuploidy data for each cancer
aneuploidy_cancer_type <- split(aneuploidy_merged_pancancer, aneuploidy_merged_pancancer$Type)

#creating counting dataframe with deletions, neutrals, and amplifications for each cancer and chromosome arm
aneuploidy_count_fast.cycling <- data.frame()
aneuploidy_count_quiescent <- data.frame()
for (i in names(aneuploidy_cancer_type)){
  for (j in 3:41){
    quiescent_deletions <- nrow(aneuploidy_merged_pancancer[which(aneuploidy_merged_pancancer$Type == i & 
                                                                    aneuploidy_merged_pancancer$Group == "Highly Quiescent" &  
                                                                    aneuploidy_merged_pancancer[,j] == -1),])
    
    quiescent_neutrals <- nrow(aneuploidy_merged_pancancer[which(aneuploidy_merged_pancancer$Type == i & 
                                                                    aneuploidy_merged_pancancer$Group == "Highly Quiescent" &  
                                                                    aneuploidy_merged_pancancer[,j] == 0),])
    
    quiescent_amplifications <- nrow(aneuploidy_merged_pancancer[which(aneuploidy_merged_pancancer$Type == i & 
                                                                    aneuploidy_merged_pancancer$Group == "Highly Quiescent" &  
                                                                    aneuploidy_merged_pancancer[,j] == +1),])
    
    fast_cycling_deletions <- nrow(aneuploidy_merged_pancancer[which(aneuploidy_merged_pancancer$Type == i & 
                                                                    aneuploidy_merged_pancancer$Group == "Fast Cycling" &  
                                                                    aneuploidy_merged_pancancer[,j] == -1),])
    
    fast_cycling_neutrals <- nrow(aneuploidy_merged_pancancer[which(aneuploidy_merged_pancancer$Type == i & 
                                                                   aneuploidy_merged_pancancer$Group == "Fast Cycling" &  
                                                                   aneuploidy_merged_pancancer[,j] == 0),])
    
    fast_cycling_amplifications <- nrow(aneuploidy_merged_pancancer[which(aneuploidy_merged_pancancer$Type == i & 
                                                                         aneuploidy_merged_pancancer$Group == "Fast Cycling" &  
                                                                         aneuploidy_merged_pancancer[,j] == +1),])
    
    quiescent_table <- cbind(Cancer= as.character(i), Chr_arm= as.character(colnames(aneuploidy_merged_pancancer)[j]), Group= "Highly Quiescent",
                             Deletions= as.numeric(quiescent_deletions), Neutral= as.numeric(quiescent_neutrals), Amplification= as.numeric(quiescent_amplifications))
    fast_cycling_table <- cbind(Cancer= as.character(i), Chr_arm= as.character(colnames(aneuploidy_merged_pancancer)[j]), Group= "Fast Cycling",
                                Deletions= as.numeric(fast_cycling_deletions), Neutral= as.numeric(fast_cycling_neutrals), Amplification= as.numeric(fast_cycling_amplifications))
    
    
    aneuploidy_count_quiescent <- rbind(aneuploidy_count_quiescent, quiescent_table)
    aneuploidy_count_fast.cycling <- rbind(aneuploidy_count_fast.cycling, fast_cycling_table)
  }
}

#merging the aneuploidy count dataframes 
aneuploidy_count <- rbind(aneuploidy_count_fast.cycling, aneuploidy_count_quiescent)
aneuploidy_count[,4:6] <- lapply(aneuploidy_count[,4:6], as.numeric)

#adding columns that show the proportion of amplifications, deletions, and neutral
aneuploidy_count <- cbind(aneuploidy_count, aneuploidy_count[, 4:6]/rowSums(aneuploidy_count[, 4:6])) 
colnames(aneuploidy_count)[7:9] <- c("Deletions_proportion", "Neutral_proportion", "Amplifications_proportion")

#plotting the proportions for each arm as stacked bars
#creating a dataframe where the proportions are long
aneuploidy_count_long <- melt(aneuploidy_count[,c(1:3, 7:9)], id.vars=c("Cancer", "Group", "Chr_arm"))

#plotting the proportion
pdf("aneuploidy_count_plots.pdf", width = 11, height = 6)
for (i in unique(aneuploidy_count_long$Cancer)){
  print(ggplot(aneuploidy_count_long[aneuploidy_count_long$Cancer == i & aneuploidy_count_long$Group == "Highly Quiescent",], 
         aes(x = factor(Chr_arm, level = unique(aneuploidy_count_long$Chr_arm)), y = value, fill = variable)) +
    geom_bar(stat = "identity") + 
    ggtitle(paste("Highly Quiescent"), i)+
    theme(legend.position="top")+
    scale_fill_manual(values=c("brown2", "grey69", "chartreuse4")))
  
  print(ggplot(aneuploidy_count_long[aneuploidy_count_long$Cancer == i & aneuploidy_count_long$Group == "Fast Cycling",], 
         aes(x = factor(Chr_arm, level = unique(aneuploidy_count_long$Chr_arm)), y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    ggtitle(paste("Fast Cycling"), i)+
    theme(legend.position="top")+
    scale_fill_manual(values=c("brown2", "grey69", "chartreuse4")))
}
dev.off()

#testing for statistical differences between groups for each arm and for each cancer
test <- cbind(aneuploidy_count_quiescent, aneuploidy_count_fast.cycling)
test[,c(4:6, 10:12)] <- lapply(test[,c(4:6, 10:12)], as.numeric)

test_dataframe_deletions <- data.frame()
test_dataframe_neutral <- data.frame()
test_dataframe_amplifications <- data.frame()
for (i in 1:nrow(test)){
  test_dataframe_deletions <- rbind(test_dataframe_deletions, c(p.value_del= fisher.test(matrix(c(test[i,4], test[i,10], sum(test[i,5]+test[i,6]), 
                                                               sum(test[i,11]+test[i,12])), ncol = 2, nrow = 2))[1],
                                                               odds.ratio_del= fisher.test(matrix(c(test[i,4], test[i,10], sum(test[i,5]+test[i,6]), 
                                                                                    sum(test[i,11]+test[i,12])), ncol = 2, nrow = 2))[[3]][[1]]))
  
  test_dataframe_neutral <- rbind(test_dataframe_neutral, c(p.value_neutral= fisher.test(matrix(c(test[i,5], test[i,11], sum(test[i,4]+test[i,6]), 
                                                                                          sum(test[i,10]+test[i,12])), ncol = 2, nrow = 2))[1],
                                                            odds.ratio_neutral= fisher.test(matrix(c(test[i,5], test[i,11], sum(test[i,4]+test[i,6]), 
                                                                                            sum(test[i,10]+test[i,12])), ncol = 2, nrow = 2))[[3]][[1]]))
  
  
  test_dataframe_amplifications <- rbind(test_dataframe_amplifications, c(p.value_ampl= fisher.test(matrix(c(test[i,6], test[i,12], sum(test[i,4]+test[i,5]), 
                                                                                   sum(test[i,10]+test[i,11])), ncol = 2, nrow = 2))[1],
                                                                          odds.ratio_ampl= fisher.test(matrix(c(test[i,6], test[i,12], sum(test[i,4]+test[i,5]), 
                                                                                                          sum(test[i,10]+test[i,11])), ncol = 2, nrow = 2))[[3]][[1]]))
}
test_dataframe <- as.data.frame(cbind(Cancer= test[,"Cancer"], Chr_arm= test[,"Chr_arm"], 
                                      adj.p.value_del= p.adjust(test_dataframe_deletions[,1], method = "BH"), Odds.ratio_del= test_dataframe_deletions[,2],
                                      adj.p.value_neutral= p.adjust(test_dataframe_neutral[,1], method = "BH"), Odds.ratio_neutral= test_dataframe_neutral[,2],
                                      adj.p.value_ampl= p.adjust(test_dataframe_amplifications[,1], method = "BH"), Odds.ratio_ampl= test_dataframe_amplifications[,2]))

test_dataframe[,3:8] <- lapply(test_dataframe[,3:8], as.numeric)

#carrying out the analysis pancancer (all cancers for each arm)
aneuploidy_count_pancancer <- cbind(cbind(aggregate(aneuploidy_count[aneuploidy_count$Group == "Highly Quiescent", 4:6], by=list(Chr_arm=aneuploidy_count[aneuploidy_count$Group == "Highly Quiescent",]$Chr_arm), FUN=sum), Group= "Highly Quiescent"),
                                    cbind(aggregate(aneuploidy_count[aneuploidy_count$Group == "Fast Cycling", 4:6], by=list(Chr_arm=aneuploidy_count[aneuploidy_count$Group == "Fast Cycling",]$Chr_arm), FUN=sum), Group= "Fast Cycling"))

test_dataframe_deletions_pancancer <- data.frame()
test_dataframe_amplifications_pancancer <- data.frame()
test_dataframe_neutral_pancancer <- data.frame()
for (i in 1:nrow(aneuploidy_count_pancancer)){
  test_dataframe_deletions_pancancer <- rbind(test_dataframe_deletions_pancancer, c(p.value_del= fisher.test(matrix(c(aneuploidy_count_pancancer[i,2], aneuploidy_count_pancancer[i,7], sum(aneuploidy_count_pancancer[i,3]+aneuploidy_count_pancancer[i,4]), 
                                                                                     sum(aneuploidy_count_pancancer[i,8]+aneuploidy_count_pancancer[i,9])), ncol = 2, nrow = 2))[1],
                                                                                     Odds.ratio_del= fisher.test(matrix(c(aneuploidy_count_pancancer[i,2], aneuploidy_count_pancancer[i,7], sum(aneuploidy_count_pancancer[i,3]+aneuploidy_count_pancancer[i,4]), 
                                                                                                          sum(aneuploidy_count_pancancer[i,8]+aneuploidy_count_pancancer[i,9])), ncol = 2, nrow = 2))[[3]][[1]]))
  
  test_dataframe_neutral_pancancer <- rbind(test_dataframe_neutral_pancancer, c(p.value_neutral= fisher.test(matrix(c(aneuploidy_count_pancancer[i,3], aneuploidy_count_pancancer[i,8], sum(aneuploidy_count_pancancer[i,2]+aneuploidy_count_pancancer[i,4]), 
                                                                                               sum(aneuploidy_count_pancancer[i,7]+aneuploidy_count_pancancer[i,9])), ncol = 2, nrow = 2))[1],
                                                                                Odds.ratio_neutral= fisher.test(matrix(c(aneuploidy_count_pancancer[i,3], aneuploidy_count_pancancer[i,8], sum(aneuploidy_count_pancancer[i,2]+aneuploidy_count_pancancer[i,4]), 
                                                                                                                               sum(aneuploidy_count_pancancer[i,7]+aneuploidy_count_pancancer[i,9])), ncol = 2, nrow = 2))[[3]][[1]]))
  
  test_dataframe_amplifications_pancancer <- rbind(test_dataframe_amplifications_pancancer, c(p.value_ampl= fisher.test(matrix(c(aneuploidy_count_pancancer[i,4], aneuploidy_count_pancancer[i,9], sum(aneuploidy_count_pancancer[i,2]+aneuploidy_count_pancancer[i,3]), 
                                                                                 sum(aneuploidy_count_pancancer[i,7]+aneuploidy_count_pancancer[i,8])), ncol = 2, nrow = 2))[1],
                                                                                 Odds.ratio_ampl= fisher.test(matrix(c(aneuploidy_count_pancancer[i,4], aneuploidy_count_pancancer[i,9], sum(aneuploidy_count_pancancer[i,2]+aneuploidy_count_pancancer[i,3]), 
                                                                                                      sum(aneuploidy_count_pancancer[i,7]+aneuploidy_count_pancancer[i,8])), ncol = 2, nrow = 2))[[3]][[1]]))
}
test_dataframe_pancancer <- cbind(test_dataframe_deletions_pancancer, test_dataframe_neutral_pancancer, test_dataframe_amplifications_pancancer)
colnames(test_dataframe_pancancer) <- c("p.value_del", "Odds_ratio_del", "p.value_neutral", "Odds_ratio_neutral", "p.value_ampl", "Odds_ratio_ampl")
test_dataframe_pancancer$adj_p.value_del <- p.adjust(as.numeric(test_dataframe_pancancer$p.value_del), method = "BH")
test_dataframe_pancancer$adj_p.value_neutral <- p.adjust(as.numeric(test_dataframe_pancancer$p.value_neutral), method = "BH")
test_dataframe_pancancer$adj_p.value_ampl <- p.adjust(as.numeric(test_dataframe_pancancer$p.value_ampl), method = "BH")
test_dataframe_pancancer$chr_arm <- aneuploidy_count_pancancer[,1]

#creating two dataframes, termed direction, where the change between groups for each cancer and each arm is recorded
#e.g. if in quiescent CESC samples chr 1p is more deleted, it will appear as -1, if more amplified +1, and if more stable 0.
direction_dataframe_quiescent <- as.data.frame(matrix(0, ncol = length(unique(test_dataframe$Chr_arm)), nrow = length(unique(test_dataframe$Cancer))))
direction_dataframe_fast.cycling <- as.data.frame(matrix(0, ncol = length(unique(test_dataframe$Chr_arm)), nrow = length(unique(test_dataframe$Cancer))))
colnames(direction_dataframe_quiescent) <- unique(test_dataframe$Chr_arm)
colnames(direction_dataframe_fast.cycling) <- unique(test_dataframe$Chr_arm)
rownames(direction_dataframe_quiescent) <- unique(test_dataframe$Cancer)
rownames(direction_dataframe_fast.cycling) <- unique(test_dataframe$Cancer)

for (i in unique(test_dataframe$Chr_arm)){
  for (j in unique(test_dataframe$Cancer)){
    if (test_dataframe[test_dataframe$Chr_arm == i & test_dataframe$Cancer == j,]$adj.p.value_del < 0.05 & test_dataframe[test_dataframe$Chr_arm == i & test_dataframe$Cancer == j,]$Odds.ratio_del > 1){
      direction_dataframe_quiescent[j, i] <- -1
    } 
    if (test_dataframe[test_dataframe$Chr_arm == i & test_dataframe$Cancer == j,]$adj.p.value_ampl < 0.05 & test_dataframe[test_dataframe$Chr_arm == i & test_dataframe$Cancer == j,]$Odds.ratio_ampl > 1){
      direction_dataframe_quiescent[j, i] <- 1
    }
  }
}

for (i in unique(test_dataframe$Chr_arm)){
  for (j in unique(test_dataframe$Cancer)){
    if (test_dataframe[test_dataframe$Chr_arm == i & test_dataframe$Cancer == j,]$adj.p.value_del < 0.05 & test_dataframe[test_dataframe$Chr_arm == i & test_dataframe$Cancer == j,]$Odds.ratio_del < 1){
      direction_dataframe_fast.cycling[j, i] <- -1
    } 
    if (test_dataframe[test_dataframe$Chr_arm == i & test_dataframe$Cancer == j,]$adj.p.value_ampl < 0.05 & test_dataframe[test_dataframe$Chr_arm == i & test_dataframe$Cancer == j,]$Odds.ratio_ampl < 1){
      direction_dataframe_fast.cycling[j, i] <- +1
    }
  }
}

#plotting a heatmap for the two direction dataframes
#removing the cancers that show no change across all chromosome arms



pdf("Heatmaps_Fast_Cycling_and_Quiescent.pdf")
par(mfrow=c(2,1))
pheatmap(as.matrix(direction_dataframe_fast.cycling), cluster_cols = FALSE, scale = "none",
         main = "A) Chromosome arm Aneuploidy in Fast Cycling cells across 24 Cancers", fontsize = 14)
pheatmap(as.matrix(direction_dataframe_quiescent), cluster_cols = FALSE, scale = "none",
         main = "B) Chromosome arm Aneuploidy in Highly Quiescent cells across 24 Cancers", fontsize = 14)
dev.off()

count_deleted_quiescent <- c()
count_amplified_quiescent <- c()
count_stable_quiescent <- c()

for (i in colnames(direction_dataframe_quiescent)){
  count_deleted_quiescent <- c(count_deleted_quiescent, length(which(direction_dataframe_quiescent[,i] == -1)))
  count_amplified_quiescent <- c(count_amplified_quiescent, length(which(direction_dataframe_quiescent[,i] == +1)))
  count_stable_quiescent <- c(count_stable_quiescent, length(which(direction_dataframe_quiescent[,i] == 0)))
}

direction_dataframe_quiescent["Times_deleted",] <- count_deleted_quiescent
direction_dataframe_quiescent["Times_amplified",] <- count_amplified_quiescent
direction_dataframe_quiescent["Times_stable",] <- count_stable_quiescent

count_deleted_fast.cycling <- c()
count_amplified_fast.cycling <- c()
count_stable_fast.cycling <- c()

for (i in colnames(direction_dataframe_fast.cycling)){
  count_deleted_fast.cycling <- c(count_deleted_fast.cycling, length(which(direction_dataframe_fast.cycling[,i] == -1)))
  count_amplified_fast.cycling <- c(count_amplified_fast.cycling, length(which(direction_dataframe_fast.cycling[,i] == +1)))
  count_stable_fast.cycling <- c(count_stable_fast.cycling, length(which(direction_dataframe_fast.cycling[,i] == 0)))
}

direction_dataframe_fast.cycling["Times_deleted",] <- count_deleted_fast.cycling
direction_dataframe_fast.cycling["Times_amplified",] <- count_amplified_fast.cycling
direction_dataframe_fast.cycling["Times_stable",] <- count_stable_fast.cycling


                              ##### ----- Linear Regression ------ #####


#converting values to factors and releveling to set reference level to 0
aneuploidy_merged_pancancer_linear_regression <- aneuploidy_merged_pancancer
aneuploidy_merged_pancancer_linear_regression[,3:41] <- lapply(aneuploidy_merged_pancancer_linear_regression[,3:41], as.factor)
aneuploidy_merged_pancancer_linear_regression[,3:41] <- lapply(aneuploidy_merged_pancancer_linear_regression[,3:41], relevel, ref= 2)

colnames(aneuploidy_merged_pancancer_linear_regression) <- c("Sample", "Type", "chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                                           "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", "chr13q", "chr14q", 
                                           "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", "chr20p", "chr20q", "chr21q", "chr22q", "X.1", 
                                           "X", "CancerType", "Barcode", "QuiescenceScore", "Group", "Regression")
aneuploidy_merged_pancancer_linear_regression <- aneuploidy_merged_pancancer_linear_regression[,-c(1:2, 42:45, 47:48)]

#building the linear model
mylinear <- lm(QuiescenceScore ~ ., data = na.omit(aneuploidy_merged_pancancer_linear_regression))

#removing insignificant terms from the model and building a new model
mylinear_only_sig <- lm(QuiescenceScore ~ chr1p + chr2p + chr2q + chr3p + chr3q + chr4q + chr6p + chr7p + chr7q + chr8p +
                          chr10p + chr10q + chr13q + chr16q + chr17p + chr17q + chr19p + chr19q + chr22q, data = na.omit(aneuploidy_merged_pancancer_linear_regression))
#comparing the reduced to the full model with anova, #test in favour of the full model, also R squared in less in reduced
anova(mylinear, mylinear_only_sig, test= "Chisq")

#improving the model by removing terms
#based on AIC, removing terms
modAIC_backward <- MASS::stepAIC(mylinear, k = 2, direction= "backward")  
#based on AIC, removing and adding terms, backward and both produce IDENTICAL models          
modAIC_default <- MASS::stepAIC(mylinear, k = 2, direction= "both")

#the AIC improved model performs better than the full model, based on ANOVA and R squared
anova(mylinear, modAIC_default, test = "Chisq")

#based on BIC, removing terms
modBIC_backward <-MASS::stepAIC(mylinear, direction = "backward", k = log(nrow(na.omit(aneuploidy_merged_pancancer_linear_regression))))
#based on BIC, removing and adding terms, backward and both produce IDENTICAL models
modBIC_default <-MASS::stepAIC(mylinear, direction = "both", k = log(nrow(na.omit(aneuploidy_merged_pancancer_linear_regression))))

#the BIC improved model performs worse than the full model, both in terms of anova and R squared
anova(mylinear, modBIC_default, test = "Chisq")
#----- Conclusion: (unbalanced) the AIC produced model performs better than the full and the BIC improved models

#balancing the datasets
aneuploidy_merged_pancancer_linear_regression_balanced <- aneuploidy_merged_pancancer
aneuploidy_merged_pancancer_linear_regression_balanced <- na.omit(aneuploidy_merged_pancancer_linear_regression_balanced)
colnames(aneuploidy_merged_pancancer_linear_regression_balanced) <- c("Sample", "Type", "chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                                                             "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", "chr13q", "chr14q", 
                                                             "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", "chr20p", "chr20q", "chr21q", "chr22q", "X.1", 
                                                             "X", "CancerType", "Barcode", "QuiescenceScore", "Group", "Regression")
#sampling random Highly Quiescent samples to match the number of Fast Cycling samples
sampled <- sample(x= aneuploidy_merged_pancancer_linear_regression_balanced[aneuploidy_merged_pancancer_linear_regression_balanced$Group == "Highly Quiescent",]$Sample, 
       size= length(aneuploidy_merged_pancancer_linear_regression_balanced[aneuploidy_merged_pancancer_linear_regression_balanced$Group == "Fast Cycling",]$Sample), 
       replace = FALSE)
#creating the sample dataframe
aneuploidy_merged_pancancer_linear_regression_balanced <- aneuploidy_merged_pancancer_linear_regression_balanced[aneuploidy_merged_pancancer_linear_regression_balanced$Sample %in% sampled | 
                                                                                                                   aneuploidy_merged_pancancer_linear_regression_balanced$Group == "Fast Cycling",]

aneuploidy_merged_pancancer_linear_regression_balanced <- aneuploidy_merged_pancancer_linear_regression_balanced[,-c(1:2, 42:45, 47:48)] 
aneuploidy_merged_pancancer_linear_regression_balanced[,1:39] <- lapply(aneuploidy_merged_pancancer_linear_regression_balanced[,1:39], as.factor)
aneuploidy_merged_pancancer_linear_regression_balanced[,1:39] <- lapply(aneuploidy_merged_pancancer_linear_regression_balanced[,1:39], relevel, ref= 2)

#building the balanced model
#the balanced model has a slightly better R squared compared to the unbalanced model
mylinear_balanced <- lm(QuiescenceScore ~ ., data = aneuploidy_merged_pancancer_linear_regression_balanced)

#removing significant terms from the balanced model
mylinear_balanced_only_sig <- lm(QuiescenceScore ~ chr2p + chr2q + chr3p + chr3q + chr7p + chr7q + chr10p + chr10q + chr13q + chr16q + chr19p + chr19q,
                                 data = aneuploidy_merged_pancancer_linear_regression_balanced)
#the full model performs better than the reduced both in terms of R2 and anova test
anova(mylinear_balanced, mylinear_balanced_only_sig, test= "Chisq")

#improving the balanced model by removing terms
#backward and default AIC functions produce IDENTICAl models
modAIC_backward_balanced <- MASS::stepAIC(mylinear_balanced, k = 2, direction= "backward")
modAIC_default_balanced <- MASS::stepAIC(mylinear_balanced, k = 2, direction= "both")

#comparing the AIC improved model to the full model, #test in favour of the reduced model, and R squared slightly improved
anova(mylinear_balanced, modAIC_default_balanced, test= "Chisq")

#based on BIC, both approaches produce IDENTICAL models
modBIC_backward_balanced <-MASS::stepAIC(mylinear_balanced, direction = "backward", k = log(nrow(na.omit(aneuploidy_merged_pancancer_linear_regression_balanced))))
modBIC_default_balanced <-MASS::stepAIC(mylinear_balanced, direction = "both", k = log(nrow(na.omit(aneuploidy_merged_pancancer_linear_regression_balanced))))

#comparing the BIC improved model to the full model, #test in favour of the full model, #R squared also reduced
anova(mylinear_balanced, modBIC_default_balanced, test= "Chisq")
#----- Conclusion: (balanced) the AIC model performs better than the full model

#----- FINAL Conclusion: (final) the balanced model performs better than the unbalanced in terms of AIC. #The AIC improved balanced model is the best.

#testing the balanced model in the samples not found in the balanced dataset
aneuploidy_merged_pancancer_linear_regression_other_half <- aneuploidy_merged_pancancer[!aneuploidy_merged_pancancer$Sample %in% sampled &
                                                                                          aneuploidy_merged_pancancer$Group != "Fast Cycling",]
aneuploidy_merged_pancancer_linear_regression_other_half <- na.omit(aneuploidy_merged_pancancer_linear_regression_other_half)

colnames(aneuploidy_merged_pancancer_linear_regression_other_half) <- c("Sample", "Type", "chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                                                                        "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", "chr13q", "chr14q", 
                                                                        "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", "chr20p", "chr20q", "chr21q", "chr22q", "X.1", 
                                                                        "X", "CancerType", "Barcode", "QuiescenceScore", "Group", "Regression")

aneuploidy_merged_pancancer_linear_regression_other_half <- aneuploidy_merged_pancancer_linear_regression_other_half[,-c(1:2, 42:45, 47:48)] 
aneuploidy_merged_pancancer_linear_regression_other_half[,1:39] <- lapply(aneuploidy_merged_pancancer_linear_regression_other_half[,1:39], as.factor)
aneuploidy_merged_pancancer_linear_regression_other_half[,1:39] <- lapply(aneuploidy_merged_pancancer_linear_regression_other_half[,1:39], relevel, ref= 2)

#predicting the quiescence score on the samples outise of the training dataset
aneuploidy_merged_pancancer_linear_regression_other_half$prediction <- predict(modAIC_default_balanced, aneuploidy_merged_pancancer_linear_regression_other_half)

#scatter plot to see how well it predicted the quiescence score
ggplot(data = aneuploidy_merged_pancancer_linear_regression_other_half, aes(x= QuiescenceScore, y = prediction)) + geom_point() +
  geom_smooth(method=lm) + ggtitle("Linear Regression balanced AIC improved model, Prediction on 900 Highly Quiescent samples") +
  stat_cor(label.x = 3, label.y = 34) +
  stat_regline_equation(label.x = 3, label.y = 32)

#testing the performance of the model based on classification
for (i in 1:nrow(aneuploidy_merged_pancancer_linear_regression_other_half)){
  ifelse(test = aneuploidy_merged_pancancer_linear_regression_other_half$QuiescenceScore[i] < 0, 
         yes = aneuploidy_merged_pancancer_linear_regression_other_half$Group_correct[i] <- "Fast Cycling",
         no = aneuploidy_merged_pancancer_linear_regression_other_half$Group_correct[i] <- "Highly Quiescent")
  
  ifelse(test = aneuploidy_merged_pancancer_linear_regression_other_half$prediction[i] < 0, 
         yes = aneuploidy_merged_pancancer_linear_regression_other_half$Group_predicted[i] <- "Fast Cycling",
         no = aneuploidy_merged_pancancer_linear_regression_other_half$Group_predicted[i] <- "Highly Quiescent")
}

#classified 86.8% of the quiescent samples correctly
nrow(aneuploidy_merged_pancancer_linear_regression_other_half[aneuploidy_merged_pancancer_linear_regression_other_half$Group_correct == 
                                                                aneuploidy_merged_pancancer_linear_regression_other_half$Group_predicted,])/
  nrow(aneuploidy_merged_pancancer_linear_regression_other_half)

#training the linear regression model in 2/3 balanced data
sampled_quiescent <- sample(x = aneuploidy_merged_pancancer_na_omit[aneuploidy_merged_pancancer_na_omit$Group == "Highly Quiescent",]$Sample, 
                            size = nrow(aneuploidy_merged_pancancer_na_omit[aneuploidy_merged_pancancer_na_omit$Group == "Fast Cycling",])*(2/3), replace = FALSE)

sampled_fast.cycling <- sample(x = aneuploidy_merged_pancancer_na_omit[aneuploidy_merged_pancancer_na_omit$Group == "Fast Cycling",]$Sample, 
                               size = nrow(aneuploidy_merged_pancancer_na_omit[aneuploidy_merged_pancancer_na_omit$Group == "Fast Cycling",])*(2/3), replace = FALSE)


aneuploidy_merged_pancancer_linear_regression_balanced_two_thirds <- aneuploidy_merged_pancancer_na_omit[aneuploidy_merged_pancancer_na_omit$Sample %in% sampled_quiescent |
                                                                                                           aneuploidy_merged_pancancer_na_omit$Sample %in% sampled_fast.cycling,]

colnames(aneuploidy_merged_pancancer_linear_regression_balanced_two_thirds) <- c("Sample", "Type", "chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                                                                        "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", "chr13q", "chr14q", 
                                                                        "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", "chr20p", "chr20q", "chr21q", "chr22q", "X.1", 
                                                                        "X", "CancerType", "Barcode", "QuiescenceScore", "Group", "Regression")

aneuploidy_merged_pancancer_linear_regression_balanced_two_thirds <- aneuploidy_merged_pancancer_linear_regression_balanced_two_thirds[,-c(1:2, 42:45, 47:48)] 
aneuploidy_merged_pancancer_linear_regression_balanced_two_thirds[,1:39] <- lapply(aneuploidy_merged_pancancer_linear_regression_balanced_two_thirds[,1:39], as.factor)
aneuploidy_merged_pancancer_linear_regression_balanced_two_thirds[,1:39] <- lapply(aneuploidy_merged_pancancer_linear_regression_balanced_two_thirds[,1:39], relevel, ref= 2)

#building the model
mylinear_balanced_trained_in_two_thirds <- lm(QuiescenceScore ~., data = aneuploidy_merged_pancancer_linear_regression_balanced_two_thirds)

#creating the testing dataset ()
aneuploidy_merged_pancancer_linear_regression_other_one_third <- aneuploidy_merged_pancancer_na_omit[!aneuploidy_merged_pancancer_na_omit$Sample %in% sampled_quiescent &
                                                 !aneuploidy_merged_pancancer_na_omit$Sample %in% sampled_fast.cycling,]

colnames(aneuploidy_merged_pancancer_linear_regression_other_one_third) <- c("Sample", "Type", "chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                                                                                 "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", "chr13q", "chr14q", 
                                                                                 "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", "chr20p", "chr20q", "chr21q", "chr22q", "X.1", 
                                                                                 "X", "CancerType", "Barcode", "QuiescenceScore", "Group", "Regression")

aneuploidy_merged_pancancer_linear_regression_other_one_third <- aneuploidy_merged_pancancer_linear_regression_other_one_third[,-c(1:2, 42:45, 47:48)] 
aneuploidy_merged_pancancer_linear_regression_other_one_third[,1:39] <- lapply(aneuploidy_merged_pancancer_linear_regression_other_one_third[,1:39], as.factor)
aneuploidy_merged_pancancer_linear_regression_other_one_third[,1:39] <- lapply(aneuploidy_merged_pancancer_linear_regression_other_one_third[,1:39], relevel, ref= 2)

#testing the model
aneuploidy_merged_pancancer_linear_regression_other_one_third$prediction <- predict(mylinear_balanced_trained_in_two_thirds, aneuploidy_merged_pancancer_linear_regression_other_one_third)

#scatter plot to see how well it predicted the quiescence score
ggplot(data = aneuploidy_merged_pancancer_linear_regression_other_one_third, aes(x= QuiescenceScore, y = prediction)) + geom_point() +
  geom_smooth(method=lm) + ggtitle("Linear Regression balanced model trained in 2/3 , Prediction on 152 FC and 1044 HQ samples") +
  stat_cor(label.x = 3, label.y = 34) +
  stat_regline_equation(label.x = 3, label.y = 32)

#testing the performance of the model based on classification
for (i in 1:nrow(aneuploidy_merged_pancancer_linear_regression_other_one_third)){
  ifelse(test = aneuploidy_merged_pancancer_linear_regression_other_one_third$QuiescenceScore[i] < 0, 
         yes = aneuploidy_merged_pancancer_linear_regression_other_one_third$Group_correct[i] <- "Fast Cycling",
         no = aneuploidy_merged_pancancer_linear_regression_other_one_third$Group_correct[i] <- "Highly Quiescent")
  
  ifelse(test = aneuploidy_merged_pancancer_linear_regression_other_one_third$prediction[i] < 0, 
         yes = aneuploidy_merged_pancancer_linear_regression_other_one_third$Group_predicted[i] <- "Fast Cycling",
         no = aneuploidy_merged_pancancer_linear_regression_other_one_third$Group_predicted[i] <- "Highly Quiescent")
}

#classified 78.3% of all (quiescent only) samples correctly
nrow(aneuploidy_merged_pancancer_linear_regression_other_one_third[aneuploidy_merged_pancancer_linear_regression_other_one_third$Group_correct == 
                                                                     aneuploidy_merged_pancancer_linear_regression_other_one_third$Group_predicted,])/
  nrow(aneuploidy_merged_pancancer_linear_regression_other_one_third)
#classified 82.3% of the quiescent samples correctly 
correct_quiescent <- aneuploidy_merged_pancancer_linear_regression_other_one_third[aneuploidy_merged_pancancer_linear_regression_other_one_third$Group_correct == "Highly Quiescent",]
correct_fast.cycling <- aneuploidy_merged_pancancer_linear_regression_other_one_third[aneuploidy_merged_pancancer_linear_regression_other_one_third$Group_correct == "Fast Cycling",]

nrow(correct_quiescent[correct_quiescent$Group_predicted == "Highly Quiescent",])/nrow(correct_quiescent)
#classified 50.7% of the fast cycling samples correctly
nrow(correct_fast.cycling[correct_fast.cycling$Group_predicted == "Fast Cycling",])/nrow(correct_fast.cycling)

#improving the model with stepAIC
modAIC_default_balanced_trained_in_two_thirds <- MASS::stepAIC(mylinear_balanced_trained_in_two_thirds, k= 2, direction = "both")

#anova in favour of new model, #R^2 in favour of new model, #AIC in favour of new model
anova(mylinear_balanced_trained_in_two_thirds, modAIC_default_balanced_trained_in_two_thirds)

#testing the new step AIC improved model
aneuploidy_merged_pancancer_linear_regression_other_one_third$prediction <- predict(mylinear_balanced_trained_in_two_thirds, aneuploidy_merged_pancancer_linear_regression_other_one_third)

#scatter plot to see how well it predicted the quiescence score
p1 <- ggplot(data = aneuploidy_merged_pancancer_linear_regression_other_one_third, aes(x= QuiescenceScore, y = prediction)) + geom_point() +
  geom_smooth(method=lm) + ggtitle("A) Linear Regression prediction on 152 FC and 1044 HQ tumours") +
  stat_cor(label.x = 3, label.y = 30, size = 6) +
  stat_regline_equation(label.x = 3, label.y = 26, size = 6) +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title=element_text(size=15, face="bold")) +
  xlab("True Quiescence Score") + ylab("Predicted Quiescence Score")

#testing the performance of the model based on classification
for (i in 1:nrow(aneuploidy_merged_pancancer_linear_regression_other_one_third)){
  ifelse(test = aneuploidy_merged_pancancer_linear_regression_other_one_third$QuiescenceScore[i] < 0, 
         yes = aneuploidy_merged_pancancer_linear_regression_other_one_third$Group_correct[i] <- "Fast Cycling",
         no = aneuploidy_merged_pancancer_linear_regression_other_one_third$Group_correct[i] <- "Highly Quiescent")
  
  ifelse(test = aneuploidy_merged_pancancer_linear_regression_other_one_third$prediction[i] < 0, 
         yes = aneuploidy_merged_pancancer_linear_regression_other_one_third$Group_predicted[i] <- "Fast Cycling",
         no = aneuploidy_merged_pancancer_linear_regression_other_one_third$Group_predicted[i] <- "Highly Quiescent")
}

#classified 76.8% of all (quiescent only) samples correctly
nrow(aneuploidy_merged_pancancer_linear_regression_other_one_third[aneuploidy_merged_pancancer_linear_regression_other_one_third$Group_correct == 
                                                                     aneuploidy_merged_pancancer_linear_regression_other_one_third$Group_predicted,])/
  nrow(aneuploidy_merged_pancancer_linear_regression_other_one_third)
#classified 81.1% of the quiescent samples correctly 
correct_quiescent <- aneuploidy_merged_pancancer_linear_regression_other_one_third[aneuploidy_merged_pancancer_linear_regression_other_one_third$Group_correct == "Highly Quiescent",]
correct_fast.cycling <- aneuploidy_merged_pancancer_linear_regression_other_one_third[aneuploidy_merged_pancancer_linear_regression_other_one_third$Group_correct == "Fast Cycling",]

nrow(correct_quiescent[correct_quiescent$Group_predicted == "Highly Quiescent",])/nrow(correct_quiescent)
#classified 47.4% of the fast cycling samples correctly
nrow(correct_fast.cycling[correct_fast.cycling$Group_predicted == "Fast Cycling",])/nrow(correct_fast.cycling)

#plot for classification
classification_summary <- rbind(data.frame(group = "Highly Quiescent", Prediction = "Correct", 
                                                     value = nrow(correct_quiescent[correct_quiescent$Group_predicted == "Highly Quiescent",]),
                                                                  percentage = 100 * nrow(correct_quiescent[correct_quiescent$Group_predicted == "Highly Quiescent",])/ nrow(correct_quiescent)), 
                                          data.frame(group = "Highly Quiescent", Prediction = "False", 
                                                     value = nrow(correct_quiescent[correct_quiescent$Group_predicted != "Highly Quiescent",]),
                                                     percentage = 100 * nrow(correct_quiescent[correct_quiescent$Group_predicted != "Highly Quiescent",])/ nrow(correct_quiescent)),
                                          data.frame(group = "Fast Cycling", Prediction = "Correct", 
                                                     value = nrow(correct_fast.cycling[correct_fast.cycling$Group_predicted == "Fast Cycling",]),
                                                     percentage = 100 * nrow(correct_fast.cycling[correct_fast.cycling$Group_predicted == "Fast Cycling",])/ nrow(correct_fast.cycling)), 
                                          data.frame(group = "Fast Cycling", Prediction = "False", 
                                                     value = nrow(correct_fast.cycling[correct_fast.cycling$Group_predicted != "Fast Cycling",]),
                                                     percentage = 100 * nrow(correct_fast.cycling[correct_fast.cycling$Group_predicted != "Fast Cycling",])/ nrow(correct_fast.cycling)))

classification_summary$percentage <- paste(round(classification_summary$percentage, digits = 1),"%")


blank_theme <- theme_minimal()+
  theme(
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    plot.title=element_text(size=15, face="bold")
  )

p2 <- ggplot(data= classification_summary[classification_summary$group == "Highly Quiescent",], aes(x="", y=value, fill= Prediction)) + 
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) + blank_theme +
  geom_label(aes(label = percentage),
             position = position_stack(vjust = 0.55),
             show.legend = FALSE, size = 6) + 
  scale_fill_manual(values= c("chartreuse4", "brown2")) + 
  ggtitle("B) Highly Quiescent tumours")

p3 <- ggplot(data= classification_summary[classification_summary$group == "Fast Cycling",], aes(x="", y=value, fill= Prediction)) + 
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) + blank_theme +
  geom_label(aes(label = percentage),
             position = position_stack(vjust = 0.55),
             show.legend = FALSE, size = 6) + 
  scale_fill_manual(values= c("chartreuse4", "brown2")) +
  ggtitle("C) Fast Cycling tumours")

grid.arrange(p1, arrangeGrob(p2, p3, ncol = 2, top = textGrob("Linear model, Group prediction", gp=gpar(fontsize= 16, font=1))), nrow = 2)




                                            ### ---------------------   LOGISTIC REGRESSION ANALYSIS   --------------------- ### 



#using logistic regression
aneuploidy_merged_pancancer_logistic_regression <- aneuploidy_merged_pancancer
aneuploidy_merged_pancancer_logistic_regression <- aneuploidy_merged_pancancer_logistic_regression[,-c(1:2, 42:47)]
aneuploidy_merged_pancancer_logistic_regression[,1:40] <- lapply(aneuploidy_merged_pancancer_logistic_regression[,1:40], as.factor)
aneuploidy_merged_pancancer_logistic_regression[,1:39] <- lapply(aneuploidy_merged_pancancer_logistic_regression[,1:39], relevel, ref = 2)

colnames(aneuploidy_merged_pancancer_logistic_regression) <- c("chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                                                    "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", "chr13q", "chr14q", 
                                                    "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", "chr20p", "chr20q", "chr21q", "chr22q", "Regression")
#building the full logistic regression model
mylogit <- glm(Regression ~ ., data= na.omit(aneuploidy_merged_pancancer_logistic_regression), family= "binomial")

#eliminating the non-significant terms from the model
mylogit_only_sig <- glm(Regression ~ chr3p + chr3q + chr4p + chr4q + chr5q + chr6p + chr7p + chr8p + chr8q + chr13q + 
                          chr16q + chr17p + chr19p + chr19q + chr20q + chr21q + chr22q, 
                        data= na.omit(aneuploidy_merged_pancancer_logistic_regression), family= "binomial")
#comparing the reduced model to the old model, #test in favour of the full model, and R squared is also better
anova(mylogit, mylogit_only_sig, test="Chisq")
pR2(mylogit)
pR2(mylogit_only_sig)

#wald test for significant terms in the model
wald_test_vector <- c()
for (i in 1:79){
  wald_test_vector <- c(wald_test_vector, wald.test(b= coef(mylogit), Sigma= vcov(mylogit), Terms = i)[[6]][[1]][3][[1]])
}
dataframe <- as.data.frame(cbind(element= rownames(as.data.frame(exp(coef(mylogit)))), Wald.test.value= as.numeric(wald_test_vector)))
rownames(dataframe) <- 1:nrow(dataframe)

#eliminating terms that failed the wald test, #model IDENTICAL to the significant terms model
mylogit_wald_test <- glm(Regression ~ chr3p + chr3q + chr4p + chr4q + chr5q + chr6p + chr7p + chr8p + chr8q + 
                           chr16q + chr17p + chr19p + chr19q + chr20q + chr21q + chr22q,
                            data= na.omit(aneuploidy_merged_pancancer_logistic_regression), family= "binomial")
#comparing the reduced model to the old model, #test in favour of the old model
anova(mylogit, mylogit_wald_test, test="Chisq")
# ---- Conclusion: (unbalanced) The full unbalanced model performs better than the reduced unbalanced models

#balancing the dataset (HQ and FC) samples 
aneuploidy_merged_pancancer_logistic_regression_na_omit <- aneuploidy_merged_pancancer
aneuploidy_merged_pancancer_logistic_regression_na_omit <- na.omit(aneuploidy_merged_pancancer_logistic_regression_na_omit)

sampled_quiescent <- sample(x = aneuploidy_merged_pancancer_na_omit[aneuploidy_merged_pancancer_na_omit$Group == "Highly Quiescent",]$Sample, 
                            size = nrow(aneuploidy_merged_pancancer_na_omit[aneuploidy_merged_pancancer_na_omit$Group == "Fast Cycling",])*(2/3), replace = FALSE)

sampled_fast.cycling <- sample(x = aneuploidy_merged_pancancer_na_omit[aneuploidy_merged_pancancer_na_omit$Group == "Fast Cycling",]$Sample, 
                               size = nrow(aneuploidy_merged_pancancer_na_omit[aneuploidy_merged_pancancer_na_omit$Group == "Fast Cycling",])*(2/3), replace = FALSE)

aneuploidy_merged_pancancer_logistic_regression_balanced <- aneuploidy_merged_pancancer_logistic_regression_na_omit[aneuploidy_merged_pancancer_logistic_regression_na_omit$Sample %in% sampled_quiescent | 
                                                                                                                      aneuploidy_merged_pancancer_logistic_regression_na_omit$Sample %in% sampled_fast.cycling,]

aneuploidy_merged_pancancer_logistic_regression_balanced <- aneuploidy_merged_pancancer_logistic_regression_balanced[,c(3:41, 48)]
aneuploidy_merged_pancancer_logistic_regression_balanced[,1:40] <- lapply(aneuploidy_merged_pancancer_logistic_regression_balanced[,1:40], as.factor)
aneuploidy_merged_pancancer_logistic_regression_balanced[,1:39] <- lapply(aneuploidy_merged_pancancer_logistic_regression_balanced[,1:39] , relevel, ref = 2)

colnames(aneuploidy_merged_pancancer_logistic_regression_balanced) <- c("chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                                                               "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", "chr13q", "chr14q", 
                                                               "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", "chr20p", "chr20q", "chr21q", "chr22q", "Regression")

#building the full balanced model,#performs better than the unbalanced based on R squared
mylogit_balanced <- glm(Regression ~ ., data= aneuploidy_merged_pancancer_logistic_regression_balanced, family= "binomial")
pR2(mylogit_balanced)

#removing insignificant terms from the balanced model
mylogit_balanced_only_sig <- glm(Regression ~ chr1p + chr1q + chr3p + chr6p + chr6q + chr8p + chr14q + chr15q +
                                    chr16q + chr17p + chr20q + chr21q + chr22q, data = na.omit(aneuploidy_merged_pancancer_logistic_regression_balanced), family= "binomial")
#evaluating the model, #test in favour of the old, full model
anova(mylogit_balanced, mylogit_balanced_only_sig, test= "Chisq")
pR2(mylogit_balanced)
pR2(mylogit_balanced_only_sig)

#removing terms based on wald test, model IDENTICAL to the isnignificant temrs model
wald_test_vector <- c()
for (i in 1:79){
  wald_test_vector <- c(wald_test_vector, wald.test(b= coef(mylogit_balanced), Sigma= vcov(mylogit_balanced), Terms = i)[[6]][[1]][3][[1]])
}
dataframe <- as.data.frame(cbind(element= rownames(as.data.frame(exp(coef(mylogit_balanced)))), Wald.test.value= as.numeric(wald_test_vector)))
dataframe$Wald.test.value <- as.numeric(dataframe$Wald.test.value)
rownames(dataframe) <- 1:nrow(dataframe)

#building the model based on wald test values, model identical to the significant terms model
mylogit_balanced_wald_test <- glm(Regression ~ chr1p + chr1q + chr3p + chr6p + chr6q + chr8p + chr14q + chr15q + 
                                    chr16q + chr17p + chr20q + chr21q + chr22q, data = na.omit(aneuploidy_merged_pancancer_logistic_regression_balanced), family= "binomial")
anova(mylogit_balanced, mylogit_balanced_wald_test, test= "Chisq")

# ----- Conclusion: (balanced) the full balanced model performs better than the reduced balanced models
# ----- Conclusion (FINAL): (final) the full balanced model performs better than the full unbalanced model, making it the best of all logistic regression models

testing <- aneuploidy_merged_pancancer_logistic_regression_na_omit[!aneuploidy_merged_pancancer_logistic_regression_na_omit$Sample %in% sampled,]
colnames(testing)[3:41] <- c("chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                                                                        "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", "chr13q", "chr14q", 
                                                                        "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", "chr20p", "chr20q", "chr21q", "chr22q")
testing[,c(3:41, 48)] <- lapply(testing[,c(3:41, 48)], as.factor)
testing[,3:41] <- lapply(testing[,3:41], relevel, ref = 2)

testing$prediction <- predict(mylogit_balanced, newdata= testing, type = "response")
testing$prediction_binary <- ifelse(testing$prediction > 0.5, yes= "Highly Quiescent", no = "Fast Cycling")

#the model classifies 73.7% of the samples correctly
nrow(testing[testing$Group == testing$prediction_binary,])/nrow(testing)

testing_quiescence <- testing[testing$Group == "Highly Quiescent",]
testing_fast.cycling <- testing[testing$Group == "Fast Cycling",]
#model classifies 80.2% of the quiescent samples correctly
nrow(testing_quiescence[testing_quiescence$Group == testing_quiescence$prediction_binary,])/
  nrow(testing_quiescence)
#model classifies 61.1% of the fast cycling samples correctly
nrow(testing_fast.cycling[testing_fast.cycling$Group == testing_fast.cycling$prediction_binary,])/
  nrow(testing_fast.cycling)


                                                  ### ---------------------   RANDOM FOREST ANALYSIS   --------------------- ###


#random forersts for analysis
aneuploidy_merged_pancancer_random_forest <- aneuploidy_merged_pancancer
aneuploidy_merged_pancancer_random_forest[,c(3:41, 47)] <- lapply(aneuploidy_merged_pancancer_random_forest[,c(3:41, 47)], as.factor)
colnames(aneuploidy_merged_pancancer_random_forest)[3:41] <- c("chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                             "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", "chr13q", "chr14q", 
                             "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", "chr20p", "chr20q", "chr21q", "chr22q")
aneuploidy_merged_pancancer_random_forest <- aneuploidy_merged_pancancer_random_forest[,c(3:41, 47)]

#finding what is the optimum mtry to use
oob.values <- c()
for (i in 1:length(3:41)){
  temp.model <- randomForest(formula= Group ~ ., data= na.omit(aneuploidy_merged_pancancer_random_forest), 
                               proximity = TRUE, type= classification, mtry = i)
  oob.values[i] <- temp.model$err.rate[nrow(temp.model$err.rate), 1]
  
}
#using (4) predictors at a time produces the least error
ggplot(data= (data.frame(values= oob.values, mtry = seq(1, 39, 1))), aes(x= mtry, y= values, label= mtry))+ 
  geom_line() + geom_point() + geom_text(hjust=-1, vjust=0) + ggtitle("unbalanced")

#building the random forest on ommited unbalanced data, #the model is sqewed towards Quiescence
randomforest <- randomForest(formula= Group ~ ., data= na.omit(aneuploidy_merged_pancancer_random_forest), 
                             proximity = TRUE, type= classification, mtry = 4, ntree= 1000)

#balancing the dataset
aneuploidy_merged_pancancer_random_forest_na_omit <- na.omit(aneuploidy_merged_pancancer)

sampled <- sample(x= aneuploidy_merged_pancancer_random_forest_na_omit[aneuploidy_merged_pancancer_random_forest_na_omit$Group == "Highly Quiescent",]$Sample, 
                  size= length(aneuploidy_merged_pancancer_random_forest_na_omit[aneuploidy_merged_pancancer_random_forest_na_omit$Group == "Fast Cycling",]$Sample), 
                  replace = FALSE)

aneuploidy_merged_pancancer_random_forest_balanced <- aneuploidy_merged_pancancer_random_forest_na_omit[aneuploidy_merged_pancancer_random_forest_na_omit$Sample %in% sampled | 
                                                                                                                aneuploidy_merged_pancancer_random_forest_na_omit$Group == "Fast Cycling",]

colnames(aneuploidy_merged_pancancer_random_forest_balanced)[3:41] <- c("chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                                                               "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", "chr13q", "chr14q", 
                                                               "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", "chr20p", "chr20q", "chr21q", "chr22q")
aneuploidy_merged_pancancer_random_forest_balanced[,c(3:41, 47)] <- lapply(aneuploidy_merged_pancancer_random_forest_balanced[,c(3:41, 47)], as.factor)
aneuploidy_merged_pancancer_random_forest_balanced <- aneuploidy_merged_pancancer_random_forest_balanced[,c(3:41, 47)]

#finding what is the optimum mtry to use
oob.values <- c()
for (i in 1:length(3:41)){
  temp.model <- randomForest(formula= Group ~ ., data= na.omit(aneuploidy_merged_pancancer_random_forest_balanced), 
                             proximity = TRUE, type= classification, mtry = i)
  oob.values[i] <- temp.model$err.rate[nrow(temp.model$err.rate), 1]
  
}
#using (9) predictors at a time produces the least error
ggplot(data= (data.frame(values= oob.values, mtry = seq(1, 39, 1))), aes(x= mtry, y= values, label= mtry))+ 
  geom_line() + geom_point() + geom_text(hjust=-0.5, vjust=0) + ggtitle("balanced")

#building the balanced random forest
randomforest_balanced <- randomForest(formula= Group ~ ., data= na.omit(aneuploidy_merged_pancancer_random_forest_balanced), 
                                                                         proximity = TRUE, type= classification, mtry= 9, ntree= 1000)

testing <- aneuploidy_merged_pancancer_random_forest_na_omit[!aneuploidy_merged_pancancer_random_forest_na_omit$Sample %in% sampled,]

colnames(testing)[3:41] <- c("chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                                                                        "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", "chr13q", "chr14q", 
                                                                        "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", "chr20p", "chr20q", "chr21q", "chr22q")
testing <- testing[,c(3:41, 47)]
testing[,1:39] <- lapply(testing[,1:39], as.factor)
testing[,1:39] <- lapply(testing[,1:39], relevel, ref= 2)
testing <- testing[!testing$Group == "Fast Cycling",]

testing$prediction <- predict(randomforest_balanced, testing)

#training the model on half of the balanced dataset samples
sampled_quiescent <- sample(x= aneuploidy_merged_pancancer_random_forest_na_omit[aneuploidy_merged_pancancer_random_forest_na_omit$Group == "Highly Quiescent",]$Sample,
                            size = 228, replace= FALSE)
sampled_fast.cycling <- sample(x= aneuploidy_merged_pancancer_random_forest_na_omit[aneuploidy_merged_pancancer_random_forest_na_omit$Group == "Fast Cycling",]$Sample,
                            size = 228, replace= FALSE)

aneuploidy_merged_pancancer_random_forest_balanced_half <- aneuploidy_merged_pancancer_random_forest_na_omit[aneuploidy_merged_pancancer_random_forest_na_omit$Sample %in% sampled_quiescent |
                                                                                                               aneuploidy_merged_pancancer_random_forest_na_omit$Sample %in% sampled_fast.cycling,]
colnames(aneuploidy_merged_pancancer_random_forest_balanced_half)[3:41] <- c("chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                                                                        "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", "chr13q", "chr14q", 
                                                                        "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", "chr20p", "chr20q", "chr21q", "chr22q")
aneuploidy_merged_pancancer_random_forest_balanced_half[,c(3:41, 47)] <- lapply(aneuploidy_merged_pancancer_random_forest_balanced_half[,c(3:41, 47)], as.factor)
aneuploidy_merged_pancancer_random_forest_balanced_half <- aneuploidy_merged_pancancer_random_forest_balanced_half[,c(3:41, 47)]

#finding what is the optimum mtry to use
oob.values <- c()
for (i in 1:length(3:41)){
  temp.model <- randomForest(formula= Group ~ ., data= na.omit(aneuploidy_merged_pancancer_random_forest_balanced), 
                             proximity = TRUE, type= classification, mtry = i)
  oob.values[i] <- temp.model$err.rate[nrow(temp.model$err.rate), 1]
  
}
#using (7) predictors at a time produces the least error
ggplot(data= (data.frame(values= oob.values, mtry = seq(1, 39, 1))), aes(x= mtry, y= values, label= mtry))+ 
  geom_line() + geom_point() + geom_text(hjust=-0.5, vjust=0)+ ggtitle("balanced, trained-in-half")

#building the model that is trained in half of the balanced dataset's samples
randomforest_balanced_trained_in_half <- randomForest(Group ~ ., data= aneuploidy_merged_pancancer_random_forest_balanced_half, 
                                                      proximity = TRUE, mtry= 7, type = classification)

#predicting the rest of the dataset's samples through this model
aneuploidy_merged_pancancer_random_forest_balanced_other_half <- aneuploidy_merged_pancancer_random_forest_na_omit[!aneuploidy_merged_pancancer_random_forest_na_omit$Sample %in% sampled_quiescent |
                                                                                                                     aneuploidy_merged_pancancer_random_forest_na_omit$Sample %in% sampled_fast.cycling,]
colnames(aneuploidy_merged_pancancer_random_forest_balanced_other_half)[3:41] <- c("chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                                                                             "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", "chr13q", "chr14q", 
                                                                             "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", "chr20p", "chr20q", "chr21q", "chr22q")
aneuploidy_merged_pancancer_random_forest_balanced_other_half[,c(3:41, 47)] <- lapply(aneuploidy_merged_pancancer_random_forest_balanced_other_half[,c(3:41, 47)], as.factor)
aneuploidy_merged_pancancer_random_forest_balanced_other_half <- aneuploidy_merged_pancancer_random_forest_balanced_other_half[,c(3:41, 47)]

aneuploidy_merged_pancancer_random_forest_balanced_other_half$prediction <- predict(randomforest_balanced_trained_in_half, aneuploidy_merged_pancancer_random_forest_balanced_other_half)

#random forest model predicts 75.6% of the samples correctly
nrow(aneuploidy_merged_pancancer_random_forest_balanced_other_half[aneuploidy_merged_pancancer_random_forest_balanced_other_half$Group == 
                                                                aneuploidy_merged_pancancer_random_forest_balanced_other_half$prediction,])/ nrow(aneuploidy_merged_pancancer_random_forest_balanced_other_half)

testing_quiescent <- aneuploidy_merged_pancancer_random_forest_balanced_other_half[aneuploidy_merged_pancancer_random_forest_balanced_other_half$Group == "Highly Quiescent",]
testing_fast.cycling <- aneuploidy_merged_pancancer_random_forest_balanced_other_half[aneuploidy_merged_pancancer_random_forest_balanced_other_half$Group == "Fast Cycling",]

#random forest model predicts 84.4% of the quiescent samples correctly
nrow(testing_quiescent[testing_quiescent$Group == testing_quiescent$prediction,])/nrow(testing_quiescent)
#random forest predicts 54.1% of the fast cycling samples correctly
nrow(testing_fast.cycling[testing_fast.cycling$Group == testing_fast.cycling$prediction,])/nrow(testing_fast.cycling)




#using Adaboost
adaboost <- boosting(Group ~., data= aneuploidy_merged_pancancer_random_forest_balanced_half, boos= TRUE)

aneuploidy_merged_pancancer_adaboost_other_half <-  aneuploidy_merged_pancancer_random_forest_balanced_other_half[,-41]
prediction <- predict(adaboost, aneuploidy_merged_pancancer_adaboost_other_half)

#plotting the number of Deletions of chr arms between groups across 24 cancers 
plot_dataframe_HQ <- rbind(chr= colnames(direction_dataframe_fast.cycling), Freq= direction_dataframe_quiescent["Times_deleted",], Group= "Highly Quiescent")
plot_dataframe_FC <- rbind(chr= colnames(direction_dataframe_fast.cycling), Freq= direction_dataframe_fast.cycling["Times_deleted",], Group= "Fast Cycling")
colnames(plot_dataframe_HQ) <- seq(1, length(colnames(direction_dataframe_fast.cycling)), 1) 
colnames(plot_dataframe_FC) <- seq(1, length(colnames(direction_dataframe_fast.cycling)), 1) 
plot_dataframe_HQ <- as.data.frame(t(plot_dataframe_HQ))
plot_dataframe_FC <- as.data.frame(t(plot_dataframe_FC))

plot_dataframe <- rbind(plot_dataframe_HQ, plot_dataframe_FC)

ggplot(data= plot_dataframe, aes(x= factor(chr, level = unique(chr)), y= as.numeric(Freq), fill= Group))+
  geom_bar(stat= "identity", position= "dodge", width= 0.7)+
  ggtitle("Deletions of chromosome arms across 24 cancers")+ ylab("Number of Cancers")+ xlab("Chromosome arm")+
  theme(legend.position="top")
  
#plotting the number of Amplifications of chr arms between groups across 24 cancers 
plot_dataframe_HQ <- rbind(chr= colnames(direction_dataframe_fast.cycling), Freq= direction_dataframe_quiescent["Times_amplified",], Group= "Highly Quiescent")
plot_dataframe_FC <- rbind(chr= colnames(direction_dataframe_fast.cycling), Freq= direction_dataframe_fast.cycling["Times_amplified",], Group= "Fast Cycling")
colnames(plot_dataframe_HQ) <- seq(1, length(colnames(direction_dataframe_fast.cycling)), 1) 
colnames(plot_dataframe_FC) <- seq(1, length(colnames(direction_dataframe_fast.cycling)), 1) 
plot_dataframe_HQ <- as.data.frame(t(plot_dataframe_HQ))
plot_dataframe_FC <- as.data.frame(t(plot_dataframe_FC))

plot_dataframe <- rbind(plot_dataframe_HQ, plot_dataframe_FC)

for (i in unique(plot_dataframe$chr)){
  if (sum(as.numeric(plot_dataframe[plot_dataframe$chr == i,]$Freq)) == 0){
    plot_dataframe <- plot_dataframe[!plot_dataframe$chr == i,]
  }
}

ggplot(data= plot_dataframe, aes(x= factor(chr, level = unique(chr)), y= as.numeric(Freq), fill= Group))+
  geom_bar(stat= "identity", position= "dodge")+
  ggtitle("Amplifications of chromosome arms across 24 cancers")+ ylab("Number of Cancers")+ xlab("Chromosome arm")+
  theme(legend.position="top")

#plotting proportions for pancancer chromosomal arms
aneuploidy_count_pancancer_long <- rbind(melt(aneuploidy_count_pancancer[1:4], id= c("Chr_arm")), 
                                         melt(aneuploidy_count_pancancer[6:10], id= c("Chr_arm")))

aneuploidy_count_pancancer_long <- rbind(aneuploidy_count_pancancer[1:5], aneuploidy_count_pancancer[6:10])
  

aneuploidy_count_pancancer_long$Deletions_proportion <- aneuploidy_count_pancancer_long[,2]/(aneuploidy_count_pancancer_long[,2]+ aneuploidy_count_pancancer_long[,3]+ aneuploidy_count_pancancer_long[,4])
aneuploidy_count_pancancer_long$Neutral_proportion <- aneuploidy_count_pancancer_long[,3]/(aneuploidy_count_pancancer_long[,2] + aneuploidy_count_pancancer_long[,3]+ aneuploidy_count_pancancer_long[,4])
aneuploidy_count_pancancer_long$Amplifications_proportion <- aneuploidy_count_pancancer_long[,4]/(aneuploidy_count_pancancer_long[,2] + aneuploidy_count_pancancer_long[,3]+ aneuploidy_count_pancancer_long[,4])

aneuploidy_count_pancancer_long<- melt(aneuploidy_count_pancancer_long[,c(1,5:8)], id.vars=c("Group", "Chr_arm"))


 p1 <- ggplot(aneuploidy_count_pancancer_long[aneuploidy_count_pancancer_long$Group == "Highly Quiescent",], 
               aes(x = factor(Chr_arm, level = unique(aneuploidy_count_long$Chr_arm)), y = value, fill = variable)) +
          geom_bar(stat = "identity") + 
          ggtitle("B) Highly Quiescent, Pancancer")+
          theme(axis.text.x = element_text(size = 15, angle = 90),
                axis.text.y = element_text(size = 15),
                legend.text = element_text(size = 16),
                title = element_text(size=15, face='bold'),
                legend.position="top",
                axis.title.x = element_text(size = 15, face = "bold"),
                axis.title.y = element_text(size = 15, face = "bold")) +
          scale_fill_manual(values=c("brown2", "grey69", "chartreuse4")) +
   xlab("Chromosome arm")
  
 p2 <- ggplot(aneuploidy_count_pancancer_long[aneuploidy_count_pancancer_long$Group == "Fast Cycling",], 
               aes(x = factor(Chr_arm, level = unique(aneuploidy_count_long$Chr_arm)), y = value, fill = variable)) +
          geom_bar(stat = "identity") + 
          ggtitle("A) Fast Cycling, Pancancer")+
          theme(axis.text.x = element_text(size = 15, angle = 90),
                axis.text.y = element_text(size = 15),
                legend.text = element_text(size = 16),
                title = element_text(size=15, face='bold'),
                legend.position="top",
                axis.title.x = element_text(size = 15, face = "bold"),
                axis.title.y = element_text(size = 15, face = "bold")) +
                
          scale_fill_manual(values=c("brown2", "grey69", "chartreuse4")) +
   xlab("Chromosome arm")


grid.arrange(p2, p1, nrow = 2)


#plotting proportions in a barplot, for pancancer
ggplot(aneuploidy_count_pancancer_long[aneuploidy_count_pancancer_long$variable == "Deletions_proportion",], aes(x= factor(Chr_arm, level = unique(aneuploidy_count_long$Chr_arm)), y= value, fill= Group))+
  geom_bar(stat= "identity",position= "dodge", width= 0.7)+
  ggtitle("Proportion of Deletions in Fast Cycling and Highly Quiescent samples, Pancancer")+ ylab("Proportion of samples")+ xlab("Chromosome arm")+
  theme(legend.position="top")

ggplot(aneuploidy_count_pancancer_long[aneuploidy_count_pancancer_long$variable == "Amplifications_proportion",], aes(x= factor(Chr_arm, level = unique(aneuploidy_count_long$Chr_arm)), y= value, fill= Group))+
  geom_bar(stat= "identity",position= "dodge", width= 0.7)+
  ggtitle("Proportion of Amplifications in Fast Cycling and Highly Quiescent samples, Pancancer")+ ylab("Proportion of samples")+ xlab("Chromosome arm")+
  theme(legend.position="top")

ggplot(aneuploidy_count_pancancer_long, mapping= aes(x= factor(Chr_arm, level = unique(aneuploidy_count_long$Chr_arm))))+
  geom_bar(data= subset(aneuploidy_count_pancancer_long, variable == "Amplifications_proportion"), mapping= aes(y = value, fill = Group), stat = "identity", position = "dodge")+
  geom_bar(data= subset(aneuploidy_count_pancancer_long, variable == "Deletions_proportion"), mapping= aes(y = -value, fill = Group), stat = "identity", position = "dodge")+
  ylim(-0.55, 0.55)+
  geom_hline(yintercept=0, linetype="dashed", size= 1)+
  ggtitle("Proportion of Amplifications and Deletions of Chromosomal arms for Fast Cycling and Highly Quiescent cells across 24 Cancers")+
  ylab("Proportion of samples")+ xlab("Chromosome arms")+
  theme(legend.position="top")

#creating a dataframe only for significant samples
test_dataframe_significant <- test_dataframe[test_dataframe$adj_p.value_del < 0.05 | test_dataframe$adj_p.value_ampl < 0.05 | test_dataframe$adj_p.value_neutral < 0.05,]

#plotting the overall number of changes
ggplot(data.frame(Number= c(153, 25, 142), CNV= c("Deletions", "Neutral", "Amplifications")), aes(x= "", y = Number, fill = CNV))+
  geom_col() +
  geom_text(aes(label = Number), position = position_stack(vjust = 0.5))+
  coord_polar(theta = "y")+
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), legend.position = "top")+
  scale_fill_manual(values=c("chartreuse4", "brown2",  "grey69"))+
  ggtitle("Number of significant changes in chromosomal arm aneuploidy in 24 cancers")
  
#counting which chromosome arm has the most significant changes between groups across cancers
chr_arm_deletions_count <- data.frame(table(test_dataframe_significant[test_dataframe_significant$adj_p.value_del < 0.05,]$chr_arm), Type= "Deletions")
chr_arm_neutral_count <- data.frame(table(test_dataframe_significant[test_dataframe_significant$adj_p.value_neutral < 0.05,]$chr_arm), Type= "Neutral")
chr_arm_amplifications_count <- data.frame(table(test_dataframe_significant[test_dataframe_significant$adj_p.value_ampl < 0.05,]$chr_arm), Type= "Amplifications")

chr_arm_count <- rbind(chr_arm_deletions_count, chr_arm_neutral_count, chr_arm_amplifications_count)
chr_arm_count <- merge(chr_arm_count, as.data.frame(aggregate(chr_arm_count$Freq, by=list(Category=chr_arm_count$Var1), FUN=sum)), by.x = "Var1", by.y = "Category")

ggplot(chr_arm_count, aes(x= factor(Var1, level = unique(aneuploidy_count$Chr_arm)), y = Freq, fill = Type))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values=c("chartreuse4", "brown2", "grey69"))+
  ggtitle("Distribution and Type of chromosomal arm aneuploidy that significantly differes between groups across 24 cancers")+
  theme(legend.position = "top")+
  geom_text(aes(label = Freq), size = 3, hjust = 0.5, vjust = 3, position = "stack")+
  geom_text(aes(label = x), size = 5, hjust = 0.5, vjust = 3)

#plotting deletions and amplification changes separately
ggplot(chr_arm_deletions_count, aes(x= factor(Var1, level = unique(aneuploidy_count_long$Chr_arm)), y= Freq))+
  geom_bar(stat = "identity")+
  ggtitle("Distribution of significant changes in chromosomal arm Deletions across 24 cancers")
  
ggplot(chr_arm_amplifications_count, aes(x= factor(Var1, level = unique(aneuploidy_count_long$Chr_arm)), y= Freq))+
  geom_bar(stat = "identity")+
  ggtitle("Distribution of significant changes in chromosomal arm Amplifications across 24 cancers")

ggplot(rbind(chr_arm_deletions_count, chr_arm_amplifications_count), aes(x= factor(Var1, level = unique(aneuploidy_count_long$Chr_arm)), y = Freq, fill = Type))+
  geom_bar(position="dodge", stat="identity")+
  scale_fill_manual(values=c("chartreuse4", "brown2"))+
  ggtitle("Distribution and Type of chromosomal arm aneuploidy that significantly differes between groups across 24 cancers")+
  theme(legend.position = "top")








