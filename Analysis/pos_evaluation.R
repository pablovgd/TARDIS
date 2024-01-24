### Evaluation of results for Metabease.
library(ggfortify)
library(mixOmics)

feat_table <- read_xlsx("S:/di04_limet_bioinformatics/PhD Pablo/Publicaties/WIP/TARDIS/data/metabease_metabolomics/positive/feature_table.csv")
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

##Adding extra analyses of correlation of areas.

auc <- read.csv2("S:/di04_limet_bioinformatics/PhD Pablo/Publicaties/WIP/TARDIS/data/metabease_metabolomics/positive/auc_table.csv", sep =",",header = FALSE)

areas_beata_1 <- readxl::read_excel("S:/di04_limet_bioinformatics/PhD Pablo/Publicaties/WIP/TARDIS/data/metabease_metabolomics/positive/areas_beata.xlsx",sheet = 1)
areas_beata_2 <- readxl::read_excel("S:/di04_limet_bioinformatics/PhD Pablo/Publicaties/WIP/TARDIS/data/metabease_metabolomics/positive/areas_beata.xlsx",sheet = 2)

areas_xcal <- cbind(areas_beata_1,areas_beata_2)

areas_xcal <- areas_xcal[,-244]

areas_xcal <- areas_xcal %>% relocate(Run, .after = Sample )

areas_xcal <- areas_xcal[,-1]

auc <- auc[,-1]

auc <- as_tibble(t(auc))

auc[1,1] = "Sample"

colnames(auc) = auc[1,]

auc <- auc[-1,]

remove_leading_zeros_and_id <- function(input_strings) {
  # Remove "ID" and leading zeros from each string
  result <- gsub("^ID|\\b0+", "", input_strings, perl = TRUE)
  
  # Return the result
  return(result)
}


colnames(areas_xcal) <- remove_leading_zeros_and_id(colnames(areas_xcal))
colnames(areas_xcal) <- remove_leading_zeros_and_id(colnames(areas_xcal))


auc[, 2:233] <- apply(auc[, 2:233], 2, as.numeric)

areas_xcal[areas_xcal == 0] <- NA

colnames(areas_xcal)[1] = "Sample"

xcal <- areas_xcal[,intersect(colnames(areas_xcal),colnames(auc))]

tardis <- auc[,intersect(colnames(areas_xcal),colnames(auc))]

result <- data.frame()

goodfeat <- read.csv2("S:/di04_limet_bioinformatics/PhD Pablo/Publicaties/WIP/TARDIS/data/metabease_metabolomics/positive/feature_table.csv",sep =",")
goodfeat <- goodfeat[which(goodfeat$rating == "Good"),]

xcal <- xcal[,intersect(colnames(xcal),goodfeat$Component)]
tardis <- tardis[,intersect(colnames(tardis),goodfeat$Component)]



for(component in 1:142){
  x <- xcal[,component]
  y <- unlist(tardis[,component])
  cortest <- cor.test(x,y)
  cor <- cortest$estimate
  p <- cortest$p.value
  testrest <- cbind(cor,p)
  result <- rbind(result,testrest)
}
