library(mzR)
library(Spectra)
library(MSnbase)
library(readxl)
library(xcms)
library(signal)
library(pracma)
library(MsExperiment)
#First create a spectra object of all files in the folder (QC's or STD for example)
fls <- list.files(path = "K:/shares/di04_limet_bioinformatics/PhD Pablo/Tapex_test_STD_Environage/pos/", pattern = ".mzML", recursive = F)
setwd("K:/shares/di04_limet_bioinformatics/PhD Pablo/Tapex_test_STD_Environage/pos/")
spectra <- Spectra(fls, source = MsBackendMzR())
spectra$centroided <- TRUE

data <- readMsExperiment(spectraFiles = fls)


#Load database etc...


masslist <- read_excel("K:/shares/di04_limet_bioinformatics/PhD Pablo/Tapex_test_STD_Environage/targlijst.xlsx")
polarity <- "positive" #you can change the polarity here, needs to be "positive" or "negative"
mode <- "metabolomics" #you can change the mode here, needs to be "lipidomics" or "metabolomics"
#With these two parameters, you can set the allowed errors for mass & retention time!
ppm <- 5 #ppm error for EIC extraction #POLAR = 5ppm / LIPIDOMICS =  10 ppm
deltaTR = 36 

masslist_negative <- masslist[grep("]-",masslist$`Ion adduct`,fixed = T),]
masslist_positive <- masslist[grep("]+",masslist$`Ion adduct`,fixed = T),]

#Remove unnecessary columns and rename

masslist_positive <- masslist_positive[,c("ID","Name","m/z-value","RT (min)")]
masslist_negative <- masslist_negative[,c("ID","Name","m/z-value","RT (min)")]
colnames(masslist_positive) <- c("ID","NAME","m/z","tr")
colnames(masslist_negative) <- c("ID","NAME","m/z","tr")

#Set RT in seconds & make numeric

masslist_positive$tr <- as.numeric(masslist_positive$tr) *60
masslist_positive$`m/z` <- as.numeric(masslist_positive$'m/z')

masslist_negative$tr <- as.numeric(masslist_negative$tr) *60
masslist_negative$`m/z` <- as.numeric(masslist_negative$'m/z')

#set ID as character

masslist_positive$ID <- as.character(masslist_positive$ID)
masslist_negative$ID <- as.character(masslist_negative$ID)

if(polarity == "positive"){
  dbData <- masslist_positive
}
if (polarity == "negative"){
  dbData <- masslist_negative
}

mzmed <- dbData$`m/z`
rtmed <- dbData$tr



mzdeltas <- sapply(mzmed, function(mzmed) mzmed*ppm/10^6)

#calculate mzrange based on delta mz
mzRanges <- cbind(as.numeric(mzmed) - mzdeltas, as.numeric(mzmed) + mzdeltas )

#if an upper m/z boundary is lower than the minimum m/z range, set it to the minimum m/z
indexTemp <- which(mzRanges[,2] < min(spectra@backend@spectraData@listData$basePeakMZ))
mzRanges[indexTemp,] <- min(spectra@backend@spectraData@listData$basePeakMZ)

#if a lower m/z boundary is higher than the max m/z, set it to the max m/z
indexTemp <- which(mzRanges[,1] > max(spectra@backend@spectraData@listData$basePeakMZ))
mzRanges[indexTemp,] <- max(spectra@backend@spectraData@listData$basePeakMZ)

#if an upper limit is higher than the max mz & the lower limit is smaller than the max m/z, set the upper limit to the max m/z
mzRanges[which(mzRanges[,2] > max(spectra@backend@spectraData@listData$basePeakMZ) & mzRanges[,1] < max(spectra@backend@spectraData@listData$basePeakMZ)),2] <- max(spectra@backend@spectraData@listData$basePeakMZ)

#if a upper limit is larger than the minimum & the lower limit is lower than the minimum, set the lower limit to the min m/z
mzRanges[which(mzRanges[,2] > min(spectra@backend@spectraData@listData$basePeakMZ) & mzRanges[,1] < min(spectra@backend@spectraData@listData$basePeakMZ)),1] <- min(spectra@backend@spectraData@listData$basePeakMZ)



rtRanges <- cbind(as.numeric(rtmed) - deltaTR/2, as.numeric(rtmed) + deltaTR/2)

#for al upper RT limits -2 lower than minimum RT, change to lower limit of min rt & upper limit to min + 10
indexTemp <- which(rtRanges[,2] - 2 < min(spectra@backend@spectraData@listData$rtime))
rtRanges[indexTemp,] <- cbind(rep(min(spectra@backend@spectraData@listData$rtime),length(indexTemp)), rep(min(spectra@backend@spectraData@listData$rtime)+10,length(indexTemp)))


#for al lower RT limits +2 higher than maximum RT, change lower limit to max -10 and upper limit to max
indexTemp <- which(rtRanges[,1] + 2 > max(spectra@backend@spectraData@listData$rtime))
rtRanges[indexTemp,] <- cbind(rep(max(spectra@backend@spectraData@listData$rtime)-10,length(indexTemp)), rep(max(spectra@backend@spectraData@listData$rtime),length(indexTemp)))


#for upper limits higher than max and lower limit lower than max, change upper limit to max
rtRanges[which(rtRanges[,2] > max(spectra@backend@spectraData@listData$rtime) & rtRanges[,1] < max(spectra@backend@spectraData@listData$rtime)),2] <- max(spectra@backend@spectraData@listData$rtime)
#for upper limits higher than min and lower limit lower than min, change lower limit to min
rtRanges[which(rtRanges[,2] > min(spectra@backend@spectraData@listData$rtime) & rtRanges[,1] < min(spectra@backend@spectraData@listData$rtime)),1] <- min(spectra@backend@spectraData@listData$rtime)


k = 2

#Filter all the files per compound
compound <- filterRt(spectra,rtRanges[k,])
compound <- filterMzRange(compound, mzRanges[k,])

#Create an object with the rt and int data for all samples for that compound
res <- spectrapply(compound, f = compound$dataOrigin, FUN = ms_data_frame, BPPARAM = SnowParam(4))

#Extract data per sample

i = 1

int <- res[[i]]$intensity
rt <- res[[i]]$rtime


#estimate the local noise

noise <- xcms:::estimateChromNoise(int, trim = 0.05,
                                   minPts = 3 )
#detect the peaks using centWave

pks <- peaksWithCentWave(int,rt,peakwidth = c(3,20),snthresh = 0, firstBaselineCheck = F, fitgauss = T, noise = noise  , integrate = 2 ,prefilter= c(1,noise),extendLengthMSW = T)

#If multiple peaks, prioritize the one closest to database RT

if(nrow(pks)>1){
  rt_database = dbData$tr[k]
  distance <- as.matrix(dist(c(rt_database,pks[,1])))[,1]
  distance <- distance[distance!=0]
  min <- which.min(as.vector(distance))
  pks <- pks[min,]
}


#Plot raw data & the found peak

plot(rt,int, type = "l") + lines(rt,int,col= "red",type="l") + 

rect(xleft = pks[2], xright = pks[ 3],
     ybottom = rep(0, 1), ytop = pks[6], col = "#0000FF7D",
     border = "#00000040")


# Detect the peak with own code
library(signal)
smoothed <- sgolayfilt(int, p = 3, n = 7)  # Adjust the polynomial order (p) and window size (n) as needed



noise <- xcms:::estimateChromNoise(smoothed, trim = 0.05,
                                   minPts = 3 )
smoothed <- smoothed - noise

smoothed[smoothed<0] <- 0

border <- find_peak_points(int)

plot(rt,smoothed) + abline(v=rt[border$left]) + abline(v = rt[border$right]) 
x= rt[border$left:border$right]
y = smoothed[border$left:border$right]
# Calculate the area under the curve
auc <- trapz(x, y)

plot(rt,smoothed, type ="l") + lines(x,y,col="red") + points(rt,int) 

mean <- (noise + min(smoothed))/2

plot(rt,smoothed) + abline(h=noise) + abline(h = min(smoothed)) + abline(h = mean)



#### TRY WITH XCMS #### 
#Extract chromatogram
bpis <- chromatogram(data, aggregationFun = "max")

#Select cmp
rtr <- rtRanges[1,]
mzr <- mzRanges[1,]


#Filter both dimensions --> no peak found
chr_raw <- chromatogram(data, mz = mzr, rt = rtr)
xchr <- findChromPeaks(chr_raw, param = CentWaveParam(snthresh = 2))

#Filter only m/z
chr_raw_2 <- chromatogram(data, mz = mzr)
plot(chr_raw)

xchr_2 <- findChromPeaks(chr_raw, param = CentWaveParam(peakwidth = c(3,20),snthresh = 3, firstBaselineCheck = T, fitgauss = T, noise = 10000  , integrate = 2 ,prefilter= c(4,10000),extendLengthMSW = T))


# TRY WITH XCMS BUT NOT ON CHROMATOGRAM DATA
raw_data <- readMSData(files = fls,
                       mode = "onDisk")

raw_data_filt <- filterMz(raw_data,mzr)

xdata <- findChromPeaks(raw_data_filt, param = CentWaveParam(peakwidth = c(3,20),snthresh = 3, firstBaselineCheck = T, fitgauss = T, noise = 10000  , integrate = 2 ,prefilter= c(4,10000),extendLengthMSW = T))

chromPeaks(xdata)
