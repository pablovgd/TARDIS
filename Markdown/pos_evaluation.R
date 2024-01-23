### Evaluation of results for Metabease.
library(ggfortify)
library(mixOmics)

feat_table <- read_xlsx("C:/Users/Administrator/Desktop/positive/feat_table.xlsx")
evaluated_targets <- read_xlsx("C:/Users/Administrator/Desktop/positive/targets.xlsx")
eval <- evaluated_targets[,c("id","rating")]
colnames(eval) <- c("ID","rating")

feat_table <- merge(feat_table,eval,by = "ID")

feat_table$AUC <- log(feat_table$AUC)
feat_table$AUC <- pareto_scale(feat_table$AUC)

feat_table$MaxInt <- log(feat_table$MaxInt)
feat_table$MaxInt <- pareto_scale(feat_table$MaxInt)


good <- feat_table[which(feat_table$rating == "Good"),]

bad <- feat_table[which(feat_table$rating == "Bad"),]

ugly <- feat_table[which(feat_table$rating == "Ambiguous"),]

summary(ugly$peak_cor)
summary(good$peak_cor)
summary(bad$peak_cor)





ggplot(feat_table, aes(x=peak_cor, y=rating)) + 
  geom_boxplot() + coord_flip()


res.aov <- aov(SNR ~ rating, data = feat_table)
# Summary of the analysis
summary(res.aov)


ggplot(feat_table, aes(x=SNR, y=rating)) + 
  geom_boxplot() + coord_flip()


res.aov <- aov(SNR ~ rating, data = feat_table)
# Summary of the analysis
summary(res.aov)

ggplot(feat_table, aes(x=AUC, y=rating)) + 
  geom_boxplot() + coord_flip()


res.aov <- aov(AUC ~ rating, data = feat_table)
# Summary of the analysis
summary(res.aov)

ggplot(feat_table, aes(x=MaxInt, y=rating)) + 
  geom_boxplot() + coord_flip()


res.aov <- aov(MaxInt ~ rating, data = feat_table)
# Summary of the analysis
summary(res.aov)


pca_res <- prcomp(feat_table[,c(3,4,5,6)], scale. = TRUE)

autoplot(pca_res, data = feat_table, colour = 'rating')

auc <- read_csv("K:/shares/di04_limet_bioinformatics/PhD Pablo/Publicaties/WIP/TARDIS/data/metabease_metabolomics/positive/auc_table.csv")


