#Dummy script to check if all my functions work with example data


source("R/functions/loadTargetedDatabase.R")
source("R/functions/massRangeCleaner.R")
source("R/functions/selectSpectra.R")
source("R/functions/findTargetedCompounds.R")

path <- "K:/shares/di04_limet_bioinformatics/PhD Pablo/Tapex_test_STD_Environage/pos/" 

database <- loadTargetedDatabase("K:/shares/di04_limet_bioinformatics/PhD Pablo/Tapex_test_STD_Environage/targlijst.xlsx","positive",c("ID","Name","m/z-value","RT (min)"))

spectra <- selectSpectra(path,"QC","my_analysis")

findTargetedCompounds(path, spectra[[1]],spectra[[2]],QC_list = list(c(1,2),c(3,4)),batch_list = list(c(1:4),c(5:8)) ,5,36,database)
