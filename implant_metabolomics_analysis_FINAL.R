library(xcms)
library(magrittr)
library(Spectra)
library(SummarizedExperiment)
library(BiocParallel)
library(stats)
library(MSnbase)
library(mzR)
library(igraph)
library(tidyverse)
library(readxl)
library(ggthemes)
library(readr)
library(stringr)
library(ggpubr)
library(ggfortify)
library(ggplot2)
library(reshape2)
library(Rcpp)
library(viridis)
library(ComplexHeatmap)
library(pROC)
library(cowplot)
library(data.table)
library(pheatmap)

#set the directory that contains the folder from github. It should include three folders: "input_files", "intermediate_files", and "figures_and_tables"
mydir <- "C:/Users/johna/Box/John Wildenthal's First Thesis Manuscript/code_for_github"

setwd(mydir)

atremoval <- read_csv("input_files/atremoval_FINAL.csv")
drainMetadata <- read_csv("input_files/drainMetadata_noPHI.csv")
xbrkey <- read_csv("input_files/xbrkey.csv")

#import negative data
rawneg <- read_csv("input_files/20240131_all_implant_samples_neg_MS2_aln_noRTalign_compounds.csv")

#positive compounds were exported in different files, but contain the same information as negative files. I apologize as this increases the complexity of the code to analyze them
rawposcompounds <- read_csv("input_files/20240116_all_implant_samples_pos_MS1_aln_rerun3_compounds.csv")
rawposgapfilling <- read_csv("input_files/20240116_all_implant_samples_pos_MS1_aln_rerun3_gap_filling.csv")
rawpospeakrating <- read_csv("input_files/20240116_all_implant_samples_pos_MS1_aln_rerun3_peak_ratings.csv")

#this is to create a unique and easy to understand "compoundID" [CID] for each compound
#this consists of the mode (pos/neg), calculated molecular weight (note: this is sometimes incorrectly automatically calculated and requires manual inspection), and retention time (minutes)
rawposcompounds["CID"] <- paste0("pos_",rawposcompounds$`Calc. MW`,"_",rawposcompounds$`RT [min]`)
rawposgapfilling["CID"] <- paste0("pos_",rawposgapfilling$`Calc. MW`,"_",rawposgapfilling$`RT [min]`)
rawpospeakrating["CID"] <- paste0("pos_",rawpospeakrating$`Calc. MW`,"_",rawpospeakrating$`RT [min]`)
rawneg["CID"] <- paste0("neg_",rawneg$`Calc. MW`,"_",rawneg$`RT [min]`)

#sanity check to ensure compounds in positive mode are in the same order. These should equal zero and thus fail an "if" check if there are any mismatches
if(sum(rawposcompounds$CID != rawposgapfilling$CID)){stop("Mismatch in positive compounds and gapfilling matrix")}
if(sum(rawposcompounds$CID != rawpospeakrating$CID)){stop("Mismatch in positive compounds and peak rating matrix")}
################################################################################
#Positive processing now involves many positive files

#get the peak areas for each compound in each sample
logareapos <- data.frame(log10(rawposcompounds[,grep("^Area:",colnames(rawposcompounds))]))
rownames(logareapos) <- rawposcompounds$CID
logareapos <- t(logareapos)
logareapos <- logareapos[rownames(logareapos)!="Area..20240116_all_implant_samples_pos_XBR090.raw..F118.",] #this sample had insufficient volume which wasn't detected until partway through the runs, we wish to use Area..20240131_all_implant_samples_pos_XBR090_rerun.raw..F350 instead (a re-run of the same sample)

#same for peak area
peakratingpos <- data.frame(rawpospeakrating[,grep("^Peak Rating",colnames(rawpospeakrating))])
rownames(peakratingpos) <- rawpospeakrating$CID
peakratingpos <- t(peakratingpos)
peakratingpos[is.na(peakratingpos)] <- 0 #one additional step to change NAs to 0 so they don't interfere with R logical operators like >, = , etc
peakratingpos <- peakratingpos[rownames(peakratingpos)!="Peak.Rating..20240116_all_implant_samples_pos_XBR090.raw..F118.",]

#same for gap fill status
gapfillstatuspos <- data.frame(rawposgapfilling[,grep("^Gap Fill Status: ",colnames(rawposgapfilling))])
rownames(gapfillstatuspos) <- rawposgapfilling$CID
gapfillstatuspos <- t(gapfillstatuspos)
gapfillstatuspos <- gapfillstatuspos[rownames(gapfillstatuspos)!="Gap.Fill.Status..20240116_all_implant_samples_pos_XBR090.raw..F118.",]

#same for gap status
gapstatuspos <- data.frame(rawposgapfilling[,grep("^Gap Status: ",colnames(rawposgapfilling))])
rownames(gapstatuspos) <- rawposgapfilling$CID
gapstatuspos <- t(gapstatuspos)
gapstatuspos <- gapstatuspos[rownames(gapstatuspos)!="Gap.Status..20240116_all_implant_samples_pos_XBR090.raw..F118.",]

################################################################################
#Negative preprocessing is the same as positive mode above

logareaneg <- data.frame(log10(rawneg[,grep("^Area:",colnames(rawneg))]))
rownames(logareaneg) <- paste0("neg_",rawneg$`Calc. MW`,"_",rawneg$`RT [min]`)
logareaneg <- t(logareaneg)
logareaneg <- logareaneg[rownames(logareaneg)!="Area..20240131_all_implant_samples_neg_XBR090.raw..F118.",] #this sample had injected by the time the positive mode data tipped us off something had happened, so we had to re-inject negative mode too and remove the original file from analysis

peakratingneg <- data.frame(rawneg[,grep("^Peak Rating",colnames(rawneg))])
rownames(peakratingneg) <- paste0("neg_",rawneg$`Calc. MW`,"_",rawneg$`RT [min]`)
peakratingneg <- t(peakratingneg)
peakratingneg[is.na(peakratingneg)] <- 0
peakratingneg <- peakratingneg[rownames(peakratingneg)!="Peak.Rating..20240131_all_implant_samples_neg_XBR090.raw..F118.",]

gapfillstatusneg <- data.frame(rawneg[,grep("^Gap Fill Status: ",colnames(rawneg))])
rownames(gapfillstatusneg) <- paste0("neg_",rawneg$`Calc. MW`,"_",rawneg$`RT [min]`)
gapfillstatusneg <- t(gapfillstatusneg)
gapfillstatusneg <- gapfillstatusneg[rownames(gapfillstatusneg)!="Gap.Fill.Status..20240131_all_implant_samples_neg_XBR090.raw..F118.",]

gapstatusneg <- data.frame(rawneg[,grep("^Gap Status: ",colnames(rawneg))])
rownames(gapstatusneg) <- paste0("neg_",rawneg$`Calc. MW`,"_",rawneg$`RT [min]`)
gapstatusneg <- t(gapstatusneg)
gapstatusneg <- gapstatusneg[rownames(gapstatusneg)!="Gap.Status..20240131_all_implant_samples_neg_XBR090.raw..F118.",]
#######################################################################################################
#restrict our analysis for now to only samples at time of implant removal
#this also orders the dataset in the same order as they are in in our metadata file
#this pulls rows from logareapos from files that have runIDs from a sample taken at removal
logareaposatremoval <- logareapos[sapply(atremoval$runID, function(i){grep(i, rownames(logareapos))}),]
peakratingposatremoval <- peakratingpos[sapply(atremoval$runID, function(i){grep(i, rownames(peakratingpos))}),]
gapfillstatusposatremoval <- gapfillstatuspos[sapply(atremoval$runID, function(i){grep(i, rownames(gapfillstatuspos))}),]
gapstatusposatremoval <- gapstatuspos[sapply(atremoval$runID, function(i){grep(i, rownames(gapstatuspos))}),]

logareanegatremoval <- logareaneg[sapply(atremoval$runID, function(i){grep(i, rownames(logareaneg))}),]
peakratingnegatremoval <- peakratingneg[sapply(atremoval$runID, function(i){grep(i, rownames(peakratingneg))}),]
gapfillstatusnegatremoval <- gapfillstatusneg[sapply(atremoval$runID, function(i){grep(i, rownames(gapfillstatusneg))}),]
gapstatusnegatremoval <- gapstatusneg[sapply(atremoval$runID, function(i){grep(i, rownames(gapstatusneg))}),]

#more sanity checks, just in case! This checks that the rownames come from the same file individually for positive and negative modes
if(sum(str_remove(rownames(logareaposatremoval),"Area..")!=str_remove(rownames(peakratingposatremoval),"Peak.Rating.."))){stop("Mismatch in positive area and peakrating matrices")}
if(sum(str_remove(rownames(logareaposatremoval),"Area..")!=str_remove(rownames(gapfillstatusposatremoval),"Gap.Fill.Status.."))){stop("Mismatch in positive compounds and gapfill status matrices")}
if(sum(str_remove(rownames(logareaposatremoval),"Area..")!=str_remove(rownames(gapstatusposatremoval),"Gap.Status.."))){stop("Mismatch in positive compounds and gap status matrices")}

if(sum(str_remove(rownames(logareanegatremoval),"Area..")!=str_remove(rownames(peakratingnegatremoval),"Peak.Rating.."))){stop("Mismatch in negative area and peakrating matrix")}
if(sum(str_remove(rownames(logareanegatremoval),"Area..")!=str_remove(rownames(gapfillstatusnegatremoval),"Gap.Fill.Status.."))){stop("Mismatch in negative gap fill status and peakrating matrix")}
if(sum(str_remove(rownames(logareanegatremoval),"Area..")!=str_remove(rownames(gapstatusnegatremoval),"Gap.Status.."))){stop("Mismatch in negative area and gap status matrix")}

#sanity check that positive and negative rows correspond to correct samples
if(sum(sapply(str_split(sapply(str_split(rownames(logareaposatremoval),"_"), function(i){i[[6]]}),".raw"), function(j){j[[1]]})!=sapply(str_split(sapply(str_split(rownames(logareanegatremoval),"_"), function(i){i[[6]]}),".raw"), function(j){j[[1]]}))){stop("Mismatch in negative and positive matrices")}

#now combine negative and positive features into one massive dataset
#since feature names have pos or neg in them, we can just use the sample rowname from either dataset. Positive has been selected for convenience
logareaall <- cbind(logareaposatremoval, logareanegatremoval)
rownames(logareaall) <- str_remove(rownames(logareaall),"Area..20240116_all_implant_samples_pos_")
rownames(logareaall) <- str_remove(rownames(logareaall),"Area..20240131_all_implant_samples_pos_")
rownames(logareaall) <- str_remove(rownames(logareaall),"_rerun") #for XBR090
rownames(logareaall) <- sapply(str_split(rownames(logareaall),".raw"), function(i){i[[1]]})

peakratingall <- cbind(peakratingposatremoval, peakratingnegatremoval)
rownames(peakratingall) <- str_remove(rownames(peakratingall),"Peak.Rating..20240116_all_implant_samples_pos_")
rownames(peakratingall) <- str_remove(rownames(logareaall),"Peak.Rating..20240131_all_implant_samples_pos_")
rownames(peakratingall) <- str_remove(rownames(peakratingall),"_rerun")
rownames(peakratingall) <- sapply(str_split(rownames(peakratingall),".raw"), function(i){i[[1]]})

gapfillstatusall <- cbind(gapfillstatusposatremoval, gapfillstatusnegatremoval)
rownames(gapfillstatusall) <- str_remove(rownames(gapfillstatusall),"Gap.Fill.Status..20240116_all_implant_samples_pos_")
rownames(gapfillstatusall) <- str_remove(rownames(gapfillstatusall),"Gap.Fill.Status..20240131_all_implant_samples_pos_")
rownames(gapfillstatusall) <- str_remove(rownames(gapfillstatusall),"_rerun")
rownames(gapfillstatusall) <- sapply(str_split(rownames(gapfillstatusall),".raw"), function(i){i[[1]]})

gapstatusall <- cbind(gapstatusposatremoval, gapstatusnegatremoval)
rownames(gapstatusall) <- str_remove(rownames(gapstatusall),"Gap.Status..20240116_all_implant_samples_pos_")
rownames(gapstatusall) <- str_remove(rownames(gapstatusall),"Gap.Status..20240131_all_implant_samples_pos_")
rownames(gapstatusall) <- str_remove(rownames(gapstatusall),"_rerun")
rownames(gapstatusall) <- sapply(str_split(rownames(gapstatusall),".raw"), function(i){i[[1]]})

if(sum(sum(rownames(logareaall)!=rownames(peakratingall)))){stop("Error in area vs peak rating")}
if(sum(sum(rownames(logareaall)!=rownames(gapfillstatusall)))){stop("Error in area vs gapfillstatus")}
if(sum(sum(rownames(logareaall)!=rownames(gapstatusall)))){stop("Error in area vs gapstatus")}

gapfillstatusall[is.na(gapfillstatusall)] <- 0 #NAs happen if there were no gaps to fill because there was always a peak (including in blank runs). Here, '0' is unknown status which actually means a peak was detected [theoretically a '1' would be 'no gap to fill' but empirically detected peaks are 1 '0' for us so it might be a glitch in Compound Discoverer software]

dim(logareaall)

#this is log10 peak areas for ALL peaks for samples taken at time of implant removal. Save this file for future use.
write.csv(logareaall, "intermediate_files/logareaall.csv") 
################################################################################
################################################################################
#this subsets the metabolomics from above, but with all fluid samples in the study
#this includes both prospective drains, and retrospective+prospective seroma fluid
#this step excludes files from patients/breasts excluded from the study
#  see Figure 1 of the paper for details; in most cases this excludes drain fluid
#  samples from breasts that had no sample collected at time of implant removal
#  and also excludes drain collections past #3 due to small n for analysis

tokeep <- xbrkey$runID

logareaposwithdrains <- logareapos[sapply(tokeep, function(i){grep(i, rownames(logareapos))}),]
peakratingposwithdrains <- peakratingpos[sapply(tokeep, function(i){grep(i, rownames(peakratingpos))}),]
gapfillstatusposwithdrains <- gapfillstatuspos[sapply(tokeep, function(i){grep(i, rownames(gapfillstatuspos))}),]
gapstatusposwithdrains <- gapstatuspos[sapply(tokeep, function(i){grep(i, rownames(gapstatuspos))}),]

logareanegwithdrains <- logareaneg[sapply(tokeep, function(i){grep(i, rownames(logareaneg))}),]
peakratingnegwithdrains <- peakratingneg[sapply(tokeep, function(i){grep(i, rownames(peakratingneg))}),]
gapfillstatusnegwithdrains <- gapfillstatusneg[sapply(tokeep, function(i){grep(i, rownames(gapfillstatusneg))}),]
gapstatusnegwithdrains <- gapstatusneg[sapply(tokeep, function(i){grep(i, rownames(gapstatusneg))}),]

#more sanity checks, just in case! This checks that the rownames come from the same file individually for positive and negative modes
if(sum(str_remove(rownames(logareaposwithdrains),"Area..")!=str_remove(rownames(peakratingposwithdrains),"Peak.Rating.."))){stop("Mismatch in positive area and peakrating matrices")}
if(sum(str_remove(rownames(logareaposwithdrains),"Area..")!=str_remove(rownames(gapfillstatusposwithdrains),"Gap.Fill.Status.."))){stop("Mismatch in positive compounds and gapfill status matrices")}
if(sum(str_remove(rownames(logareaposwithdrains),"Area..")!=str_remove(rownames(gapstatusposwithdrains),"Gap.Status.."))){stop("Mismatch in positive compounds and gap status matrices")}

if(sum(str_remove(rownames(logareanegwithdrains),"Area..")!=str_remove(rownames(peakratingnegwithdrains),"Peak.Rating.."))){stop("Mismatch in negative area and peakrating matrix")}
if(sum(str_remove(rownames(logareanegwithdrains),"Area..")!=str_remove(rownames(gapfillstatusnegwithdrains),"Gap.Fill.Status.."))){stop("Mismatch in negative gap fill status and peakrating matrix")}
if(sum(str_remove(rownames(logareanegwithdrains),"Area..")!=str_remove(rownames(gapstatusnegwithdrains),"Gap.Status.."))){stop("Mismatch in negative area and gap status matrix")}

#sanity check that positive and negative rows correspond to correct samples
if(sum(sapply(str_split(sapply(str_split(rownames(logareaposwithdrains),"_"), function(i){i[[6]]}),".raw"), function(j){j[[1]]})!=sapply(str_split(sapply(str_split(rownames(logareanegwithdrains),"_"), function(i){i[[6]]}),".raw"), function(j){j[[1]]}))){stop("Mismatch in negative and positive matrices")}

logareawithdrains <- cbind(logareaposwithdrains, logareanegwithdrains)
rownames(logareawithdrains) <- str_remove(rownames(logareawithdrains),"Area..20240116_all_implant_samples_pos_")
rownames(logareawithdrains) <- str_remove(rownames(logareawithdrains),"Area..20240131_all_implant_samples_pos_")
rownames(logareawithdrains) <- str_remove(rownames(logareawithdrains),"_rerun") #for XBR090
rownames(logareawithdrains) <- sapply(str_split(rownames(logareawithdrains),".raw"), function(i){i[[1]]})

peakratingwithdrains <- cbind(peakratingposwithdrains, peakratingnegwithdrains)
rownames(peakratingwithdrains) <- str_remove(rownames(peakratingwithdrains),"Peak.Rating..20240116_all_implant_samples_pos_")
rownames(peakratingwithdrains) <- str_remove(rownames(logareawithdrains),"Peak.Rating..20240131_all_implant_samples_pos_")
rownames(peakratingwithdrains) <- str_remove(rownames(peakratingwithdrains),"_rerun")
rownames(peakratingwithdrains) <- sapply(str_split(rownames(peakratingwithdrains),".raw"), function(i){i[[1]]})

gapfillstatuswithdrains <- cbind(gapfillstatusposwithdrains, gapfillstatusnegwithdrains)
rownames(gapfillstatuswithdrains) <- str_remove(rownames(gapfillstatuswithdrains),"Gap.Fill.Status..20240116_all_implant_samples_pos_")
rownames(gapfillstatuswithdrains) <- str_remove(rownames(gapfillstatuswithdrains),"Gap.Fill.Status..20240131_all_implant_samples_pos_")
rownames(gapfillstatuswithdrains) <- str_remove(rownames(gapfillstatuswithdrains),"_rerun")
rownames(gapfillstatuswithdrains) <- sapply(str_split(rownames(gapfillstatuswithdrains),".raw"), function(i){i[[1]]})

gapstatuswithdrains <- cbind(gapstatusposwithdrains, gapstatusnegwithdrains)
rownames(gapstatuswithdrains) <- str_remove(rownames(gapstatuswithdrains),"Gap.Status..20240116_all_implant_samples_pos_")
rownames(gapstatuswithdrains) <- str_remove(rownames(gapstatuswithdrains),"Gap.Status..20240131_all_implant_samples_pos_")
rownames(gapstatuswithdrains) <- str_remove(rownames(gapstatuswithdrains),"_rerun")
rownames(gapstatuswithdrains) <- sapply(str_split(rownames(gapstatuswithdrains),".raw"), function(i){i[[1]]})

if(sum(sum(rownames(logareawithdrains)!=rownames(peakratingwithdrains)))){stop("Error in area vs peak rating")}
if(sum(sum(rownames(logareawithdrains)!=rownames(gapfillstatuswithdrains)))){stop("Error in area vs gapfillstatus")}
if(sum(sum(rownames(logareawithdrains)!=rownames(gapstatuswithdrains)))){stop("Error in area vs gapstatus")}

gapfillstatuswithdrains[is.na(gapfillstatuswithdrains)] <- 0 #NAs happen if there were no gaps to fill because there was always a peak (including in blank runs). Here, '0' is unknown status which actually means a peak was detected [theoretically a '1' would be 'no gap to fill' but empirically detected peaks are 1 '0' for us so it might be a glitch in Compound Discoverer software]

dim(logareawithdrains)

#this is log10 peak areas for ALL peaks for drain fluid AND seroma fluid. Save these files.
write.csv(logareawithdrains, "intermediate_files/logareawithdrains.csv") 

################################################################################
#this is the code that was used to manually integrate AD1 and AD3 levels
#these were not detected by our initial peak calls
#the mzML files are too large/numerous to fit on GitHub, but are available by request
#this code is just to show you what steps we took in integrating our data

#note that while the nomenclature in the paper was standardized to HNP1/2/3,
#  the code I generated used the nomenclature AD1/2/3 as the term "alpha defensins"
#  is more common in the surgical literature

#list.files("ALLMZML", pattern="*.mzML") 
#filelist <- paste0("ALLMZML/",list.files("ALLMZML", pattern="*.mzML"))
#alldata <- readMSData(files = filelist, mode="onDisk", verbose=TRUE)
#saveRDS(alldata, "alldata_ondisk.rds") #I advise saving this, as reading in LCMS data can take a while

alldata <- readRDS("C:/Users/johna/Box/Research Notes/IDX data/20240123 all implant samples/alldata_ondisk.rds")
allpos <- filterFile(alldata, grep("samples_pos",fileNames(alldata)))
allneg <- filterFile(alldata, grep("samples_neg",fileNames(alldata)))

instudy <- data.frame(logareawithdrains)

xbrtoget <- rownames(instudy)
xbrtoget[xbrtoget=="XBR090"] <- "XBR090_rerun" #the first injection had an error, so we need to ensure we get the re-run for this

whichfiles <- sapply(xbrtoget, function(i){grep(i, fileNames(allpos))})

#filter files for the RT/mz window of interest
subFiles <- filterFile(allpos, unlist(whichfiles))
ad1_4plus_subFiles <- subFiles
ad1_4plus_subFiles <- filterRt(ad1_4plus_subFiles, rt=c(610,630))
ad1_4plus_subFiles <- filterMz(ad1_4plus_subFiles, mz=c(861.3859-5/1000000*861.3859, 861.3859+5/1000000*861.3859))

#manually integrate AD1 abundance using manualChromPeaks
manpeakmat1 <- matrix(c(861.3859-5/1000000*861.3859, 861.3859+5/1000000*861.3859, 610,630), nrow=1)
colnames(manpeakmat1) <- c("mzmin","mzmax","rtmin","rtmax")
manAD1 <- manualChromPeaks(ad1_4plus_subFiles, chromPeaks=manpeakmat1)


mymat1 <- chromPeaks(manAD1)
vals1 <- data.frame(mymat1[,c("into","sample")])
vals1["filename"] <- sapply(vals1["sample"], function(i){fileNames(ad1_4plus_subFiles)[i]})
#this text removal may be different based on file location. Basically, it is turning the file name into the runID
vals1["xbr"] <- str_remove(vals1$filename, fixed('C:\\Users\\johna\\Box\\Research Notes\\IDX data\\20240123 all implant samples\\ALLMZML\\20240116_all_implant_samples_pos_'))
vals1$xbr <- str_remove(vals1$xbr, ".mzML")

#log-transform values and place them into the correct row of "instudy"
#Additionally, interpolate zero values to be 0.001 less than the minimum detected value of all samples (i.e. estimated limit of detection). 
#This is to avoid zero values in our dataset, which complicates certain forms of analysis (such as log ratios)
instudy["AD1"] <- sapply(rownames(instudy), function(i){ifelse(is.element(i,vals1$xbr),log10(vals1$into[vals1$xbr==i]),log10(min(vals1$into))-0.001)}) 


#####AD3
#Do the same for AD3 as was done above for AD1. Note difference in mz and RT
ad3_4plus_subFiles <- subFiles
ad3_4plus_subFiles <- filterRt(ad3_4plus_subFiles, rt=c(620,640))
ad3_4plus_subFiles <- filterMz(ad3_4plus_subFiles, mz=c(872.3826-5/1000000*872.3826, 872.3826+5/1000000*872.3826))

manpeakmat3 <- matrix(c(872.3826-5/1000000*872.3826, 872.3826+5/1000000*872.3826, 620,640), nrow=1)
colnames(manpeakmat3) <- c("mzmin","mzmax","rtmin","rtmax")

manAD3 <- manualChromPeaks(ad3_4plus_subFiles, chromPeaks=manpeakmat3)


mymat3 <- chromPeaks(manAD3)
vals3 <- data.frame(mymat3[,c("into","sample")])
vals3["filename"] <- sapply(vals3["sample"], function(i){fileNames(ad3_4plus_subFiles)[i]})
vals3["xbr"] <- str_remove(vals3$filename, fixed('C:\\Users\\johna\\Box\\Research Notes\\IDX data\\20240123 all implant samples\\ALLMZML\\20240116_all_implant_samples_pos_'))
vals3$xbr <- str_remove(vals3$xbr, ".mzML")

#log-transform values and place them into the correct row of "instudy"
#Additionally, interpolate zero values to be 0.001 less than the minimum detected value of all samples (i.e. estimated limit of detection). 
#This is to avoid zero values in our dataset, which complicates certain forms of analysis (such as log ratios)
instudy["AD3"] <- sapply(rownames(instudy), function(i){ifelse(is.element(i,vals3$xbr),log10(vals3$into[vals3$xbr==i]),log10(min(vals3$into))-0.001)})

#this file now has AD1 and AD3 manually integrated
write.csv(instudy, "intermediate_files/logareawithdrains_withADs_FINAL.csv")

################################################################################
################################################################################
##########################################################################################################################################
#for our initial analysis, we decided to focus on high-quality features present in at least 10 samples
#"high-quality" here means a peak rating >=5 in at least 10/82 samples taken at time of implant removal, as calculated by CompoundDiscoverer (Thermo)
#we assumed that features not meeting this distinction would not be good diagnositc markers
PR5IN10 <- colSums(peakratingall >= 5) >= 10 #find features with peak rating >=5 in at least 10 samples

#subset all data matrices for only those peaks
logareaPR5IN10 <- logareaall[,PR5IN10]
peakratingPR5IN10 <- peakratingall[,PR5IN10]
gapfillstatusPR5IN10 <- gapfillstatusall[,PR5IN10]
gapstatusPR5IN10 <- gapstatusall[,PR5IN10]

#this creates a list of three vectors, whether the feature was detected in positive/negative mode, the calculated MW, and the RT
getPNMWRT <- lapply(1:3, function(i){sapply(str_split(colnames(peakratingPR5IN10), "_"), function(qw){qw[[i]]})})
myposneg <- getPNMWRT[[1]]
mymws <- as.numeric(getPNMWRT[[2]])
myrts <- as.numeric(getPNMWRT[[3]])

#this is from Thermo documentation and is copied here for convenience
levs <- c(0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)
labs <- c("Unknown status",
          "No gap to fill",
          "Unable to fill",
          "Filled by arbitrary value",
          "Filled by trace area",
          "Filled by simulated peak",
          "Filled by spectrum noise",
          "Filled by matching ion",
          "Filled by re-detected peak",
          "Imputed by low area value",
          "Imputed by group median",
          "Imputed by Random Forest",
          "Skipped")

# 0 is Unknown Status which is not gap filled, so positive for compound
# 8 is Filled by Trace Area which is ~probably positive, or at least good enough
# 32 is Filled by Spectrum noise which is negative for the compound
# 64 is Filled by matching ion which is ~probably positive for a compound
# 128 is Filled by re-detected peak which is positive

#this is a binary representation of whether there is any peak/noise present or not. "Nothing present" here is the same as "Filled by Spectrum Noise" from CompoundDiscoverer
gfsYNPR5IN10 <- gapfillstatusPR5IN10!=32


#this is a correlation matrix of all features to all other features, but it calculates correlations using pairs of samples where neither has the gap status "Filled by Spectrum Noise"
newcormat <- matrix(-999, nrow=dim(logareaPR5IN10)[2], ncol=dim(logareaPR5IN10)[2])
newtmp <- lapply(1:dim(logareaPR5IN10)[2], function(qw){print(qw);sapply(1:dim(logareaPR5IN10)[2], function(qe){cor(logareaPR5IN10[gfsYNPR5IN10[,qw]&gfsYNPR5IN10[,qe],qw],logareaPR5IN10[gfsYNPR5IN10[,qw]&gfsYNPR5IN10[,qe],qe])})})

for(i in 1:dim(logareaPR5IN10)[2]){
  newcormat[i,] <- newtmp[[i]]
}
newcormat[is.na(newcormat)] <- 0

#this finds the median per-sample difference in intensity between the two peaks. It isn't used for further calculations but was used to get an idea of which peaks were larger
medIMinusJ <- matrix(-999, nrow=dim(logareaPR5IN10)[2], ncol=dim(logareaPR5IN10)[2])
tmpiminusj <- lapply(1:dim(logareaPR5IN10)[2], function(qw){print(qw);sapply(1:dim(logareaPR5IN10)[2], function(qe){median(logareaPR5IN10[,qw]-logareaPR5IN10[,qe])})})

for(i in 1:dim(logareaPR5IN10)[2]){
  medIMinusJ[i,] <- tmpiminusj[[i]]
}

#this gets the RT difference between peaks
rtdiffs <- matrix(-999, nrow=dim(logareaPR5IN10)[2], ncol=dim(logareaPR5IN10)[2])
tmprts <- lapply(1:dim(logareaPR5IN10)[2], function(qw){print(qw);abs(myrts[qw]-myrts)})

for(i in 1:dim(logareaPR5IN10)[2]){
  rtdiffs[i,] <- tmprts[[i]]
}

#this makes a network where nodes are features and edges are features that have a "Conditional correlation" >=0.9 and are within 0.1 minutes of each other
#we suspect these to be source decay families or errors in adduct calling, and remove them to avoid what are essentially duplicate rows for further statistical analysis
isf_family <- (newcormat >= 0.9) & (rtdiffs <= 0.1)
colnames(isf_family) <- colnames(logareaPR5IN10)
rownames(isf_family) <- colnames(logareaPR5IN10)

#this sets up an igraph network based on the above matrix, with metadata about correlation, median intensity, and RT difference
mygraphdf <- data.frame(which(isf_family,arr.ind = TRUE))
colnames(mygraphdf) <- c("row","col")
mygraphdf["node1"] <- sapply(mygraphdf$row, function(i){colnames(isf_family)[i]})
mygraphdf["node2"] <- sapply(mygraphdf$col, function(i){colnames(isf_family)[i]})
mygraphdf["newcormat"] <- as.numeric(sapply(1:dim(mygraphdf)[1], function(i){newcormat[mygraphdf$row[i],mygraphdf$col[i]]}))
mygraphdf["medIminusJ"] <- sapply(1:dim(mygraphdf)[1], function(i){medIMinusJ[mygraphdf$row[i],mygraphdf$col[i]]})
mygraphdf["rtdiffs"] <- sapply(1:dim(mygraphdf)[1], function(i){rtdiffs[mygraphdf$row[i],mygraphdf$col[i]]})
mygraphdf <- mygraphdf[mygraphdf$row!=mygraphdf$col,]

#this sets up the associate node metadata for the above network
nodedf <- data.frame(cbind(1:dim(logareaPR5IN10)[2],colnames(logareaPR5IN10)))
colnames(nodedf) <- c("rowCol","label")
nodedf["mode"] <- sapply(str_split(nodedf$label,"_"), function(i){i[[1]]})
nodedf["MW"] <- as.numeric(sapply(str_split(nodedf$label,"_"), function(i){i[[2]]}))
nodedf["rt"] <- as.numeric(sapply(str_split(nodedf$label,"_"), function(i){i[[3]]}))
nodedf["med_intensity"] <- colMedians(logareaPR5IN10)
nodedf["mean_intensity"] <- colMeans(logareaPR5IN10)

nodedf["mz"] <- NULL
nodedf$mz[nodedf$mode=="pos"] <- sapply(nodedf$label[nodedf$mode=="pos"], function(i){rawposcompounds$`m/z`[rawposcompounds$CID==i]})
nodedf$mz[nodedf$mode=="neg"] <- sapply(nodedf$label[nodedf$mode=="neg"], function(i){rawneg$`m/z`[rawneg$CID==i]})

#this actually creates the network
mygraph <- graph_from_data_frame(mygraphdf, directed=FALSE, vertices=nodedf)

mycomponents <- components(mygraph)

#this shows that there are an estimated 2540 features after removing suspected source decay families and missed adducts/isotopes
mycomponents$no #2540 compounds


membershipdf <- data.frame(cbind(colnames(logareaPR5IN10), mycomponents$membership, colMeans(logareaPR5IN10)))
colnames(membershipdf) <- c("name","cluster","meanValue")
membershipdf$cluster <- as.numeric(membershipdf$cluster)
membershipdf$meanValue <- as.numeric(membershipdf$meanValue)
membershipdf[order(membershipdf$cluster),]

#this finds the sample with the highest mean intensity of a suspected source decay family, and returns that
representativesample <- sapply(unique(membershipdf$cluster), function(i){clusterdf <- membershipdf[membershipdf$cluster==i,];
return(clusterdf$name[clusterdf$meanValue==max(clusterdf$meanValue)])})


outputPR5IN10 <- logareaPR5IN10[,representativesample]

regoutput <- outputPR5IN10
rownames(regoutput)[rownames(regoutput)=="XBR090_rerun"] <- "XBR090"

regoutput <- data.frame(cbind(rownames(regoutput), #runID of the sample
                              unlist(sapply(rownames(regoutput), function(i){atremoval$manualYN[atremoval$runID==i]})), #infection status
                              regoutput))

colnames(regoutput)[1:2] <- c("samplename","Status")
rownames(regoutput) <- NULL

withNIWCoutput <- regoutput


withNIWCoutput$Status[withNIWCoutput$Status=="NIWC"] <- "N" #for our purposes, we treat non-infectious removals as "uninfected"

#this is the dataset that was used for heatmap/volcano plot/SelEnergyPermR based analysis
write_csv(data.frame(withNIWCoutput),"intermediate_files/FINALDATA_WITH_NIWC_UNFILTERED.csv") 


################################################################################
#This is the analysis with SelEnergyPermR to get DCV scores for features that have high diagnostic value in a compositional setting

inputdata <- read_csv("intermediate_files/FINALDATA_WITH_NIWC_UNFILTERED.csv")
gettopnodes <- function(myfilename, topn=60, dcvinput=inputdata){
  print(myfilename)
  dcv_full <- readRDS(myfilename)
  dcvScores = dcv_full$lrs
  rm(dcv_full)#relatively big so dump from memory when you can
  gc()
  # drop all negative scores, those are actively unhelpful ratios
  dcvScores$rowmean[dcvScores$rowmean<0] = 0
  # scale postive scores to range [0,1]
  trainData.md = caret::preProcess(data.frame(dcvScores$rowmean),
                                   method = "range",rangeBounds = c(0,1) )
  scaledScore = stats::predict(trainData.md, data.frame(dcvScores$rowmean))
  # create list of edges, one for each ratio
  el = data.frame(Ratio = dcvScores$Ratio,Score = scaledScore[,1])
  el["node1"] <- sapply(str_split(el$Ratio,"___"), function(i){i[[1]]})
  el["node2"] <- sapply(str_split(el$Ratio,"___"), function(i){i[[2]]})

  # build logratio network
  g = igraph::graph_from_edgelist(as.matrix(el[,2:3]))
  igraph::E(g)$weight = el$Score
  # save network to file (optional)
  #igraph::write_graph(g, "output/dcv_logratio_network.graphml", format = "graphml") 
  # calculate node strength
  g_str <- igraph::strength(g)
  igraph::V(g)$strength <- g_str
  nodeStrength <- data.frame(Node = names(g_str), Str = g_str) %>%
    dplyr::arrange(dplyr::desc(Str))
  top_nodes <- nodeStrength[1:topn, ]
  top_nodes["file"] <- myfilename
  top_nodes["rank"] <- 1:topn
  return(top_nodes)
}

#these files take ~70-100+GB of memory to generate. See the attached scripts for an example of how to generate these on a SLURM cluster
#here, I have already generated a spreadsheet of the top nodes from each of 25 runs with different random seeds, and we will use that
#I have the commented code below to show how those were generated from the RData output files from the SLURM cluster we used

#alltopnodes <- lapply(paste0("UNFILTERED_DICOVAR/output/final_mdl_dcv_scores_UNFILTERED_DICOVAR_iteration_",c(1:25),".RData"), gettopnodes)
#nodedf <- do.call(rbind, alltopnodes)
#rownames(nodedf) <- NULL
#write_csv(nodedf,"topnodes_df_allmodels_unfilteredDCV.csv")
nodedf <- data.frame(read_csv("input_files/topnodes_df_allmodels_unfilteredDCV.csv")) #this is the result from the intermediate SelEnergyPerm output files


nodedfrank20 <- nodedf[nodedf$rank<=20,] #get top 20 ranked nodes in each randomized run
sort(table(nodedfrank20$Node),decreasing=TRUE)[1:20] #see which nodes are most common
sortednode20tbl <- sort(table(nodedfrank20$Node),decreasing=TRUE)
ntop20df <- data.frame(cbind(names(sortednode20tbl), sortednode20tbl)) #put nodes into a table along with # times selected out of 25 randomized runs
rownames(ntop20df) <- NULL
colnames(ntop20df) <- c("feature","times_selected")


finaldata <- read_csv("intermediate_files/FINALDATA_WITH_NIWC_UNFILTERED.csv")
#this is useful metadata for the output table
aucList <- lapply(3:dim(finaldata)[2], function(i){roc(finaldata$Status, unlist(finaldata[,i]),direction="<", levels=c("N","Y"))})
aucVal <- sapply(aucList, function(i){i$auc}) #AUC
aucLower <- sapply(aucList, function(i){ci(i)[1]}) #lower bound AUC
aucUpper <- sapply(aucList, function(i){ci(i)[3]}) #upper bound AUC
mannwhitney <- lapply(3:dim(finaldata)[2], function(i){wilcox.test(unlist(finaldata[finaldata$Status=="N",i]), unlist(finaldata[finaldata$Status=="Y",i]))}) #Wilcoxon Rank Sum test, aka Mann Whitney U Test
pMW <- sapply(mannwhitney, function(i){i$p.value}) # get just the P value
fdrMW <- p.adjust(pMW, method="BH") #adjust P values using BH method for FDR correction

aucdf <- data.frame(cbind(colnames(finaldata)[3:length(colnames(finaldata))], aucVal, aucLower, aucUpper,pMW,fdrMW))
colnames(aucdf)[1] <- "cid"

#make sure all columns are numeric so you can do addition/multiplication/etc
aucdf <- as.data.frame(aucdf)
aucdf$aucVal <- as.numeric(aucdf$aucVal)
aucdf$aucLower <- as.numeric(aucdf$aucLower)
aucdf$aucUpper <- as.numeric(aucdf$aucUpper)
aucdf$pMW <- as.numeric(aucdf$pMW)
aucdf$fdrMW <- as.numeric(aucdf$fdrMW)

#get the mean difference of feature intensity between infected and uninfected samples (mean intensity in infected samples minus mean intensity in uninfected samples)
aucdf["mean_diff"] <- sapply(3:dim(finaldata)[2], function(i){mean(unlist(finaldata[finaldata$Status=="Y",i]))-mean(unlist(finaldata[finaldata$Status=="N",i]))})
aucdf["log2_mean_diff"] <- aucdf$mean_diff / log10(2) #change log10 differences to log2 differences, which are slightly more human-readable
aucdf["fold_change"] <- 10^aucdf$mean_diff #the fold-change in infected vs uninfected samples, which is even more human-readable

aucdf["neglog10FDR"] <- -log10(aucdf$fdrMW)

#get associated metadata (AUC, P, FDR, etc) for each of the top nodes we selected using DCV scores and output them into a table
top20nodesoutput <- ntop20df
top20nodesoutput["AUC"] <- round(sapply(top20nodesoutput$feature, function(i){aucdf$aucVal[aucdf$cid==i]}),2)
top20nodesoutput["p (Wilcoxon)"] <- signif(sapply(top20nodesoutput$feature, function(i){aucdf$pMW[aucdf$cid==i]}),2)
top20nodesoutput["FDR (BH)"] <- signif(sapply(top20nodesoutput$feature, function(i){aucdf$fdrMW[aucdf$cid==i]}),2)
top20nodesoutput["log2 Mean Difference"] <- signif(sapply(top20nodesoutput$feature, function(i){aucdf$log2_mean_diff[aucdf$cid==i]}),2)
top20nodesoutput["Fold Change"] <- round(sapply(top20nodesoutput$feature, function(i){aucdf$fold_change[aucdf$cid==i]}),1)
top20nodesoutput["Percentage feature in top 20 rank"] <- paste0(100*(as.numeric(top20nodesoutput$times_selected)/25),"%")

top20nodesoutput <- top20nodesoutput[as.numeric(top20nodesoutput$times_selected) >=12.5,]

top20nodesoutput["Calc. MW"] <- sapply(top20nodesoutput$feature, function(i){ifelse(substr(i,1,3)=="pos",
                                                                                    rawposcompounds$`Calc. MW`[rawposcompounds$CID==i],
                                                                                    rawneg$`Calc. MW`[rawposcompounds$CID==i])})

top20nodesoutput["m/z"] <- sapply(top20nodesoutput$feature, function(i){ifelse(substr(i,1,3)=="pos",
                                                                               rawposcompounds$`m/z`[rawposcompounds$CID==i],
                                                                               rawneg$`m/z`[rawposcompounds$CID==i])})
top20nodesoutput["RT"] <- sapply(top20nodesoutput$feature, function(i){ifelse(substr(i,1,3)=="pos",
                                                                              rawposcompounds$`RT [min]`[rawposcompounds$CID==i],
                                                                              rawneg$`RT [min]`[rawposcompounds$CID==i])})

t20toprint <- top20nodesoutput[,c("feature","Calc. MW","m/z","RT","Percentage feature in top 20 rank","AUC","p (Wilcoxon)","FDR (BH)","log2 Mean Difference","Fold Change")]
write_csv(t20toprint, "intermediate_files/unfiltered_selenergyperm_top20_table.csv") #this ends up as table 2

##########################################################################################################################################

#these are used in the dotplots to estimate the limit of detection based on the mean value of five blanks that were aligned with our samples
posblanks <- data.frame(logareapos[grepl("MeOH", rownames(logareapos)),])
negblanks <- data.frame(logareaneg[grepl("MeOH",rownames(logareaneg)),])

write.csv(posblanks, "intermediate_files/posblanks.csv") 
write.csv(negblanks, "intermediate_files/negblanks.csv") 

##########################################################################################################################################


######################################################################################################################
######################################################################################################################

######################################################################################################################
######################################################################################################################
#MAKE FINAL FIGURES NOW THAT INTERMEDIATE FILES ARE MADE
#YOU SHOULD BE ABLE TO RUN THIS CODE INDEPENDENTLY OF THE ABOVE ONCE ALL INTERMEDIATE FILES ARE GENERATED

xbrkey <- read_csv("input_files/xbrkey.csv")
posblanks <- data.frame(fread("intermediate_files/posblanks.csv"), row.names=1)
negblanks <- data.frame(fread("intermediate_files/negblanks.csv"), row.names=1)
logareaall <- data.frame(fread("intermediate_files/logareaall.csv"), row.names=1)
rownames(logareaall)[rownames(logareaall)=="XBR090_rerun"] <- "XBR090"
atremoval <- read_csv("input_files/atremoval_FINAL.csv")
logareaall["XBR"] <- rownames(logareaall)

#these are metadata columns that may be useful for different plots down below
logareaall["overallYN"] <- sapply(logareaall$XBR, function(i){atremoval$overallInfStatus[atremoval$runID==i]})
logareaall["YNwithBreakdown"] <- sapply(logareaall$XBR, function(i){ifelse(atremoval$manualYN[atremoval$runID==i]!="N",atremoval$manualYN[atremoval$runID==i],ifelse(atremoval$contralateral_to_inf_sameremovaldate[atremoval$runID==i],"CONTRA","RTX"))})

logareaall$overallYN <- factor(logareaall$overallYN, levels=c("N","Y"))
logareaall$YNwithBreakdown <- factor(logareaall$YNwithBreakdown, levels = c("RTX","CONTRA","NIWC","Y"))
logareaall$InfUninf <- factor(sapply(logareaall$overallYN, function(i){ifelse(i=="Y","Infected","Uninfected")}),levels=c("Uninfected","Infected"))

#this has all samples (drain and seroma) and includes manually integrated AD1 and AD3
instudy <- data.frame(fread("intermediate_files/logareawithdrains_withADs_FINAL.csv"),row.names=1)
#this file is only samples taken at time of implant removal. We can manually add their AD1 and AD3 values from the above spreadsheet for future dotplots/box plots
logareaallADs <- logareaall 
logareaallADs["AD1"] <- sapply(rownames(logareaallADs), function(i){instudy$AD1[rownames(instudy)==i]})
logareaallADs["AD3"] <- sapply(rownames(logareaallADs), function(i){instudy$AD3[rownames(instudy)==i]})
logareaallADs["AD2"] <- logareaall$pos_1685.23651_10.418 #for convenience of plotting. Do NOT perform statistics on this duplicate column
posblanks["AD1"] <- min(instudy$AD1) #this is the estimated limit of detection calculated above based on the minimum peak area in any sample
posblanks["AD3"] <- min(instudy$AD3) #this is the estimated limit of detection calculated above based on the minimum peak area in any sample

logareaall["InfNotInf"] <- logareaall$overallYN
logareaall$InfNotInf

#graphing stuff
mytheme <- theme_classic()+theme(text = element_text(family = "Helvetica", size = 8),
                                 plot.title = element_text(family = "Helvetica", size = 8),
                                 plot.subtitle = element_text(family = "Helvetica", size = 8),
                                 axis.title = element_text(family = "Helvetica", size = 8),
                                 legend.title = element_text(family = "Helvetica", size = 8),
                                 legend.text = element_text(family = "Helvetica", size = 8),
                                 plot.caption = element_text(family = "Helvetica", size = 8),
                                 axis.text.x = element_text(family = "Helvetica", size = 8),
                                 axis.text.y = element_text(family = "Helvetica", size = 8)
)
windowsFonts("Helvetica" = windowsFont("Helvetica"))

#this will be used for the bulk of figures. It makes dotplots and boxplots of samples in infected vs uninfected breasts at time of implant removal
monochromeDotplotYN <- function(CID, mydf=logareaall, mytitle=CID,myPval="notCalculated", myblanks=posblanks){
  minval <- min(mydf[,CID])
  maxval <- max(mydf[,CID])
  rangeCID <- abs(diff(range(mydf[,CID])))
  
  p2<-ggplot(mydf, aes(x=InfUninf, y=.data[[CID]]), shape=19) + 
    mytheme+
    geom_boxplot(fill=NA, color=c("blue","red"), outliers=FALSE, linewidth=0.7)+
    geom_dotplot(binaxis='y', stackdir='center',binwidth=1/40*abs(diff(range(mydf[,CID]))))+
    ggtitle(mytitle)+ylab("log10 peak area")+theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank())+xlab("")+
    scale_y_continuous(limits = c(minval-0.2*rangeCID, maxval+0.2*rangeCID), breaks = seq(floor(min(mydf[,CID])-0.2*rangeCID),ceiling(max(mydf[,CID])+0.2*rangeCID), by = 1))+
    #theme(aspect.ratio = 2/1, plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))+
    geom_hline(yintercept=mean(myblanks[,CID]), linetype='dashed', col = "black")+
    geom_signif(
      comparisons = list(c("Uninfected", "Infected")),
      map_signif_level = TRUE,
      color="black",
      annotation=myPval,
      size = 0.5,
      vjust=-0.1,
      textsize=2.8224,
      family="Helvetica") # Add pairwise comparisons p-value
  
  plot(p2)
  return(p2)
}

#this provides annotations for the level of significant. p>0.05 NS; 0.01<=p<0.05 *; 0.001<=p<0.01 **; p<0.001 ***
plotsignifstars <- function(p){
  output <- ifelse(p>0.05,"NS",ifelse(p>0.01,"*",ifelse(p>0.001,"**","***")))
  return(output)
}
############################################################################################################
#FIGURE 5

drainMetadata <- read_csv("input_files/drainMetadata_noPHI.csv")
meltedDM <- reshape2::melt(drainMetadata[,c("runID","sampleStatus","days_before_swelling","days_before_erythema","days_before_diagnosis","days_before_removal")],id.vars = c("runID","sampleStatus"))


#make horizontal dotplot/boxplot for each drain collection (1,2,3) that progressed to infection ("Y") using ggplot
f5a <- ggplot(meltedDM[meltedDM$sampleStatus=="drain_1_Y",], aes(x=value, y=variable)) + 
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "grey", linewidth=0.7)+
  geom_boxplot(fill=NA, color="grey40", outliers=FALSE, linewidth=0.7)+
  geom_dotplot(binaxis='x', stackdir='center', binwidth=1)+
  xlab("Days since collection #1 (Day 0)")+
  ylab("")+
  mytheme+
  scale_x_continuous(limits = c(-25, 90), breaks=seq(-15,90,15))+
  scale_y_discrete(labels=c("days_before_swelling"="Swelling (n=11)","days_before_erythema"="Erythema (n=9)","days_before_diagnosis"="Diagnosis (n=16)","days_before_removal"="Implant Removal (n=16)"))

f5b <- ggplot(meltedDM[meltedDM$sampleStatus=="drain_2_Y",], aes(x=value, y=variable)) + 
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "grey", linewidth=0.7)+
  geom_boxplot(fill=NA, color="grey40", outliers=FALSE, linewidth=0.7)+
  geom_dotplot(binaxis='x', stackdir='center', binwidth=1)+
  xlab("Days since collection #2 (Day 0)")+
  ylab("")+
  mytheme+
  scale_x_continuous(limits = c(-25, 90), breaks=seq(-15,90,15))+
  scale_y_discrete(labels=c("days_before_swelling"="Swelling (n=11)","days_before_erythema"="Erythema (n=9)","days_before_diagnosis"="Diagnosis (n=15)","days_before_removal"="Implant Removal (n=15)"))
f5c <- ggplot(meltedDM[meltedDM$sampleStatus=="drain_3_Y",], aes(x=value, y=variable)) + 
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "grey", linewidth=0.7)+
  geom_boxplot(fill=NA, color="grey40", outliers=FALSE, linewidth=0.7)+
  geom_dotplot(binaxis='x', stackdir='center', binwidth=1)+
  xlab("Days since collection #3 (Day 0)")+
  ylab("")+
  mytheme+
  scale_x_continuous(limits = c(-25, 90), breaks=seq(-15,90,15))+
  scale_y_discrete(labels=c("days_before_swelling"="Swelling (n=6)","days_before_erythema"="Erythema (n=5)","days_before_diagnosis"="Diagnosis (n=9)","days_before_removal"="Implant Removal (n=9)"))

f5title <- ggdraw() + draw_label("Timing of clinical events in relation to drain fluid collections", size = 8, fontfamily = "Helvetica")

library(cowplot)
tiff("figures_and_tables/figure_5_FINAL.tiff", width=4,height=3.3, units="in",res=600)
ggdraw() +
  draw_plot(f5title, x=0.03, y=0.9, width=1, height=0.1)+
  draw_plot(f5a, x = 0.03, y = 0.6, width = 1, height = 0.3) +
  draw_plot(f5b, x = 0.03, y = 0.3, width = 1, height = 0.3) +
  draw_plot(f5c, x = 0.03, y = 0, width = 1, height = 0.3) +
  draw_plot_label(label = c("","A", "B", "C"), size = 15,
                  x = c(0,0, 0, 0), y = c(1,0.9, 0.6, 0.3))
dev.off()

sdv_5a <- drainMetadata[drainMetadata$sampleStatus=="drain_1_Y",c("runID","sampleStatus","days_before_swelling","days_before_erythema","days_before_diagnosis","days_before_removal")]
sdv_5b <- drainMetadata[drainMetadata$sampleStatus=="drain_2_Y",c("runID","sampleStatus","days_before_swelling","days_before_erythema","days_before_diagnosis","days_before_removal")]
sdv_5c <- drainMetadata[drainMetadata$sampleStatus=="drain_3_Y",c("runID","sampleStatus","days_before_swelling","days_before_erythema","days_before_diagnosis","days_before_removal")]

write.csv(sdv_5a, "supporting_data_values/sdv_5a.csv")
write.csv(sdv_5b, "supporting_data_values/sdv_5b.csv")
write.csv(sdv_5c, "supporting_data_values/sdv_5c.csv")
#put files in format for supporting data values file

############################################################################################################
#FIGURE 2
finaldata <- read_csv("intermediate_files/FINALDATA_WITH_NIWC_UNFILTERED.csv")

#generate AUCs, P Values, FDRs similarly to the above, but this time for use in the heatmap/volcano plots of Figure 3
aucList <- lapply(3:dim(finaldata)[2], function(i){roc(finaldata$Status, unlist(finaldata[,i]),direction="<", levels=c("N","Y"))})
aucVal <- sapply(aucList, function(i){i$auc})
aucLower <- sapply(aucList, function(i){ci(i)[1]})
aucUpper <- sapply(aucList, function(i){ci(i)[3]})
mannwhitney <- lapply(3:dim(finaldata)[2], function(i){wilcox.test(unlist(finaldata[finaldata$Status=="N",i]), unlist(finaldata[finaldata$Status=="Y",i]))})
pMW <- sapply(mannwhitney, function(i){i$p.value})
fdrMW <- p.adjust(pMW, method="BH")

aucdf <- data.frame(cbind(colnames(finaldata)[3:length(colnames(finaldata))], aucVal, aucLower, aucUpper,pMW,fdrMW))
colnames(aucdf)[1] <- "cid"

aucdf <- as.data.frame(aucdf)
aucdf$aucVal <- as.numeric(aucdf$aucVal)
aucdf$aucLower <- as.numeric(aucdf$aucLower)
aucdf$aucUpper <- as.numeric(aucdf$aucUpper)
aucdf$pMW <- as.numeric(aucdf$pMW)
aucdf$fdrMW <- as.numeric(aucdf$fdrMW)
aucdf["mean_diff"] <- sapply(3:dim(finaldata)[2], function(i){mean(unlist(finaldata[finaldata$Status=="Y",i]))-mean(unlist(finaldata[finaldata$Status=="N",i]))})
aucdf["log2_mean_diff"] <- aucdf$mean_diff / log10(2)

aucdf["neglog10FDR"] <- -log10(aucdf$fdrMW)

sum(aucdf$fdrMW < 0.01 & aucdf$log2_mean_diff>=2) #115 features positively correlated with infection (4-fold increase, FDR < 0.01)
sum(aucdf$fdrMW < 0.01 & aucdf$log2_mean_diff<=-2) #30 features negatively correlated with infection (4-fold decrease, FDR < 0.01)

poscorCIDs <- aucdf$cid[aucdf$fdrMW < 0.01 & aucdf$log2_mean_diff>=2] #get the compound IDs for positive correlates
negcorCIDs <- aucdf$cid[aucdf$fdrMW < 0.01 & aucdf$log2_mean_diff<=-2] #get the compound IDs for negative correlates

#make the volcano plot of log2 mean difference vs -log10(FDR)
myvolcano <- ggplot(aucdf, aes(x=log2_mean_diff, y=neglog10FDR)) + 
  geom_vline(xintercept = c(2, -2), linetype = "dashed", color = "grey")+
  geom_hline(yintercept = c(2), linetype = "dashed", color = "grey")+
  geom_point(size=0.7) +
  mytheme +
  xlab("mean log2-fold change in infection")+
  ylab("-log10 FDR-adjusted p value")+
  scale_x_continuous(breaks = seq(-4,4,2), limits=c(-5,5)) +
  scale_y_continuous(breaks = seq(0,8,2))

library(ComplexHeatmap)
col = list(`Infection Status` = c("Uninfected" = "#56B4E9", "Infected" = "#E69F00"))

ha <- HeatmapAnnotation(
  `Infection Status` = ifelse(finaldata$Status=="Y","Infected","Uninfected"),
  col = col,
  show_annotation_name = FALSE,
  annotation_legend_param = list(`Infection Status` = list( title_gp = gpar(fontsize = 8), # Change font size for annotation legend title 
                                                            labels_gp = gpar(fontsize = 8))) # Change font size for annotation legend labels
)

scaleddata <- t(scale(as.matrix(finaldata[,c(-1,-2)])))
#make the heatmap using ComplexHeatmap package
ht <- Heatmap(scaleddata, 
              name = "Z-score", #title of legend
              column_title = "Sample", row_title = "Feature",
              show_row_names = FALSE, show_column_names = FALSE,
              use_raster=FALSE, top_annotation = ha,
              column_title_gp = gpar(fontsize = 8), # Change font size for column names 
              row_title_gp = gpar(fontsize = 8),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 8), # Change font size for legend title 
                                          labels_gp = gpar(fontsize = 8))
)

grobht <- grid.grabExpr(draw(ht))

tiff("figures_and_tables/figure_2_FINAL.tiff", width = 5, height = 5, units = "in", res = 600)
ggdraw() +
  draw_plot(grobht, x = 0, y = 1/2, width = 1, height = 1/2) +
  draw_plot(myvolcano, x = 0, y = 0, width = 0.75, height = 1/2) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0), y = c(1, 1/2))

dev.off()

############################################################################################################
#FIGURE 3


#make dotplots for all identified compounds selected by DCV scores from intermediate_files/unfiltered_selenergyperm_top20_table.csv
diacspDot <- monochromeDotplotYN("pos_143.11832_1.314", myPval = plotsignifstars(aucdf$fdrMW[aucdf$cid=="pos_143.11832_1.314"]), mytitle = "Diacetylspermine")
acspDot <- monochromeDotplotYN("pos_244.22611_1.022", myPval = plotsignifstars(aucdf$fdrMW[aucdf$cid=="pos_244.22611_1.022"]), mytitle = "Acetylspermine")
hexsphDot <- monochromeDotplotYN("pos_461.33469_13.882", myPval = plotsignifstars(aucdf$fdrMW[aucdf$cid=="pos_461.33469_13.882"]), mytitle = "Hexosyl-sphingosine")
TrpGluDot <- monochromeDotplotYN("pos_333.13215_6.697", myPval = plotsignifstars(aucdf$fdrMW[aucdf$cid=="pos_333.13215_6.697"]), mytitle = "Trp-Glu")
SerLeuDot <- monochromeDotplotYN("pos_218.12658_4.867", myPval = plotsignifstars(aucdf$fdrMW[aucdf$cid=="pos_218.12658_4.867"]), mytitle = "Ser-Leu")
ValValDot <- monochromeDotplotYN("pos_216.14734_4.603", myPval = plotsignifstars(aucdf$fdrMW[aucdf$cid=="pos_216.14734_4.603"]), mytitle = "Val-Val")
PL1Dot <- monochromeDotplotYN("pos_556.32152_5.537", myPval = plotsignifstars(aucdf$fdrMW[aucdf$cid=="pos_556.32152_5.537"]), mytitle = "Peptide-like\n556.3 Da")
PL2Dot <- monochromeDotplotYN("pos_1077.58075_5.559", myPval = plotsignifstars(aucdf$fdrMW[aucdf$cid=="pos_1077.58075_5.559"]), mytitle = "Peptide-like\n1077.6 Da")



tiff(filename = "figures_and_tables/figure_3_FINAL.tiff", width = 3, height = 6, units = "in", res = 600)
ggdraw() +
  draw_plot(diacspDot, x = 0, y = 3/4, width = 0.5, height = 1/4) +
  draw_plot(acspDot, x = 0.5, y = 3/4, width = 0.5, height = 1/4) +
  draw_plot(hexsphDot, x = 0, y = 1/2, width = 0.5, height = 1/4) +
  draw_plot(TrpGluDot, x = 0.5, y = 1/2, width = 0.5, height = 1/4) +
  draw_plot(SerLeuDot, x = 0, y = 1/4, width = 0.5, height = 1/4) +
  draw_plot(ValValDot, x = 0.5, y = 1/4, width = 0.5, height = 1/4) +
  draw_plot(PL1Dot, x = 0, y = 0, width = 0.5, height = 1/4) +
  draw_plot(PL2Dot, x = 0.5, y = 0, width = 0.5, height = 1/4) +
  draw_plot_label(label = c("A", "B", "C","D","E","F","G","H"), size = 15,
                  x = c(0, 0.5, 0, 0.5, 0, 0.5, 0, 0.5), y = c(1, 1, 3/4, 3/4, 1/2, 1/2, 1/4, 1/4))
dev.off()

sdv_3a <- logareaall[,c("InfUninf","pos_143.11832_1.314")]
sdv_3b <- logareaall[,c("InfUninf","pos_244.22611_1.022")]
sdv_3c <- logareaall[,c("InfUninf","pos_461.33469_13.882")]
sdv_3d <- logareaall[,c("InfUninf","pos_333.13215_6.697")]
sdv_3e <- logareaall[,c("InfUninf","pos_218.12658_4.867")]
sdv_3f <- logareaall[,c("InfUninf","pos_216.14734_4.603")]
sdv_3g <- logareaall[,c("InfUninf","pos_556.32152_5.537")]
sdv_3h <- logareaall[,c("InfUninf","pos_1077.58075_5.559")]

write.csv(sdv_3a, "supporting_data_values/sdv_3a.csv")
write.csv(sdv_3b, "supporting_data_values/sdv_3b.csv")
write.csv(sdv_3c, "supporting_data_values/sdv_3c.csv")
write.csv(sdv_3d, "supporting_data_values/sdv_3d.csv")
write.csv(sdv_3e, "supporting_data_values/sdv_3e.csv")
write.csv(sdv_3f, "supporting_data_values/sdv_3f.csv")
write.csv(sdv_3g, "supporting_data_values/sdv_3g.csv")
write.csv(sdv_3h, "supporting_data_values/sdv_3h.csv")
############################################################################################################
#FIGURE 4

#compound IDs for Figure 4 for features that were manually searched for (AD2, clTyr, brTyr) or manually integrated (AD1, AD3)
AD2CID <- "pos_1685.23651_10.418"
cltyrCID <- "pos_215.03473_5.386"
brTyrCID <- "pos_258.98419_5.868"
AD1CID <- "AD1"
AD3CID <- "AD3"

#get the correlation of AD1 vs AD3 and AD1 vs AD2
ad1vs2cor <- signif(cor(logareaallADs[,AD1CID],logareaallADs[,AD2CID]),2)
ad1vs3cor <- signif(cor(logareaallADs[,AD1CID],logareaallADs[,AD3CID]),2)

#note that p values here are unadjusted as they weren't part of the initial 2540 feature screening above and were manually searched for
AD1dot <- monochromeDotplotYN(AD1CID,mytitle="HNP1",mydf=logareaallADs, myPval = plotsignifstars(wilcox.test(unlist(logareaallADs[logareaallADs$InfUninf=="Infected",AD1CID]),unlist(logareaallADs[logareaallADs$InfUninf=="Uninfected",AD1CID]))$p.value))
AD2dot <- monochromeDotplotYN(AD2CID,mytitle="HNP2",mydf=logareaallADs, myPval = plotsignifstars(wilcox.test(unlist(logareaallADs[logareaallADs$InfUninf=="Infected",AD2CID]),unlist(logareaallADs[logareaallADs$InfUninf=="Uninfected",AD2CID]))$p.value))
AD3dot <- monochromeDotplotYN(AD3CID,mytitle="HNP3",mydf=logareaallADs, myPval = plotsignifstars(wilcox.test(unlist(logareaallADs[logareaallADs$InfUninf=="Infected",AD3CID]),unlist(logareaallADs[logareaallADs$InfUninf=="Uninfected",AD3CID]))$p.value))
cltyrdot <- monochromeDotplotYN(cltyrCID,mydf=logareaallADs, mytitle="3-chlorotyrosine", myPval = plotsignifstars(wilcox.test(unlist(logareaallADs[logareaallADs$InfUninf=="Infected",cltyrCID]),unlist(logareaallADs[logareaallADs$InfUninf=="Uninfected",cltyrCID]))$p.value))

#this makes correlation plots so you can see AD1 and AD2 have high correlation, with AD3 having more deviations
AD1vsAD2 <- ggplot(logareaallADs, aes(x=AD1, y=AD2)) + geom_point(size=0.7) +mytheme +xlab("HNP1")+ylab("HNP2")+ggtitle(paste0("HNP1 vs HNP2 (R=",ad1vs2cor,")"))+theme(plot.title = element_text(hjust = 0.5))
AD1vsAD3 <- ggplot(logareaallADs, aes(x=AD1, y=AD3)) + geom_point(size=0.7) +mytheme +xlab("HNP1")+ylab("HNP3")+ggtitle(paste0("HNP1 vs HNP3 (R=",ad1vs3cor,")"))+theme(plot.title = element_text(hjust = 0.5))

#these are the AUC values reported in the paper, but don't fit into any specific table
roc(logareaallADs$InfUninf, logareaallADs$AD1) #AUC 0.93
roc(logareaallADs$InfUninf, logareaallADs$AD2) #AUC 0.92
roc(logareaallADs$InfUninf, logareaallADs$AD3) #AUC 0.78
roc(logareaallADs$InfUninf, logareaallADs$pos_215.03473_5.386) #AUC 0.72
roc(logareaallADs$InfUninf, logareaallADs$pos_258.98419_5.868) #AUC 0.68



tiff(filename = "figures_and_tables/figure_4_FINAL.tiff", width = 3, height = 6*3/4, units = "in", res = 600)
ggdraw() +
  draw_plot(AD1dot, x = 0, y = 2/3, width = 0.5, height = 1/3) +
  draw_plot(AD2dot, x = 0.5, y = 2/3, width = 0.5, height = 1/3) +
  draw_plot(AD3dot, x = 0, y = 1/3, width = 0.5, height = 1/3) +
  draw_plot(cltyrdot, x = 0.5, y = 1/3, width = 0.5, height = 1/3) +
  draw_plot(AD1vsAD2, x = 0, y = 0, width = 0.5, height = 1/3) +
  draw_plot(AD1vsAD3, x = 0.5, y = 0, width = 0.5, height = 1/3) +
  draw_plot_label(label = c("A", "B", "C","D","E","F"), size = 15,
                  x = c(0, 0.5, 0,0.5,0,0.5), y = c(1, 1, 2/3,2/3,1/3,1/3))
dev.off()

sdv_4a <- logareaallADs[,c("InfUninf","AD1")]
sdv_4b <- logareaallADs[,c("InfUninf","pos_1685.23651_10.418")]
sdv_4c <- logareaallADs[,c("InfUninf","AD3")]
sdv_4d <- logareaallADs[,c("InfUninf","pos_215.03473_5.386")]


write.csv(sdv_4a, "supporting_data_values/sdv_4a.csv")
write.csv(sdv_4b, "supporting_data_values/sdv_4b.csv")
write.csv(sdv_4c, "supporting_data_values/sdv_4c.csv")
write.csv(sdv_4d, "supporting_data_values/sdv_4d.csv")


############################################################################################################
#FIGURE 6 AND TABLE 4
atremoval["patient_side"] <- paste0(atremoval$patient,"_",atremoval$side)

xbrkey["patient_side"] <- paste0(xbrkey$patient,"_",xbrkey$side)
xbrkey["drainNum"] <- paste0(xbrkey$sampletype,"_",xbrkey$collectionNum)
xbrkey["drainStatus"] <- xbrkey$drainNum

#sanity check to verify all plotted drains are from the first three collections, not truly needed
xbrkey["finalInfStatus"] <- sapply(xbrkey$patient_side, function(i){ifelse(is.element(i,atremoval$patient_side),atremoval$overallInfStatus[atremoval$patient_side==i],"NOTINSTUDY")})
xbrkey["finalDrainStatus"] <- paste0(xbrkey$drainStatus,"_",xbrkey$finalInfStatus)
xbrkey["includeForLongitudinal"] <- is.element(xbrkey$runID, atremoval$runID) | is.element(xbrkey$finalDrainStatus, c("drain_1_N","drain_1_Y","drain_2_N","drain_2_Y","drain_3_N","drain_3_Y")) 

instudyLongitudinal <- instudy #create a separate version that we can alter
instudyLongitudinal <- instudyLongitudinal[sapply(rownames(instudyLongitudinal),function(i){xbrkey$includeForLongitudinal[xbrkey$runID==i]}),]
instudyLongitudinal["sampleStatus"] <- "NONE"
#for samples at time of implant removal (in atremoval data object), replace the sample status with the final infection status
instudyLongitudinal$sampleStatus[is.element(rownames(instudyLongitudinal),atremoval$runID)] <- sapply(rownames(instudyLongitudinal)[is.element(rownames(instudyLongitudinal),atremoval$runID)], function(i){atremoval$overallInfStatus[atremoval$runID==i]})
#for remaining (drain) samples, replace sample status with "drain_collectionNumber_infectionStatus" i.e. drain_1_Y is drain collection #1 from breasts that go on to have TEs removed due to infection
instudyLongitudinal$sampleStatus[instudyLongitudinal$sampleStatus=="NONE"] <- sapply(rownames(instudyLongitudinal)[instudyLongitudinal$sampleStatus=="NONE"], function(i){xbrkey$finalDrainStatus[xbrkey$runID==i]})

instudyLongitudinal$sampleStatus <- factor(instudyLongitudinal$sampleStatus, levels=c("drain_1_N","drain_1_Y","drain_2_N","drain_2_Y","drain_3_N","drain_3_Y","N","Y"))


##### MAKING TABLE 4 #####
#function to add mean and median values of a feature
addmeanmed <- function(mysampletype, mychart){
  mychart[paste0(mysampletype,"_med")] <- sapply(mychart$molecule, function(i){median(instudyLongitudinal[instudyLongitudinal$sampleStatus==mysampletype,i])})
  mychart[paste0(mysampletype,"_mean")] <- sapply(mychart$molecule, function(i){mean(instudyLongitudinal[instudyLongitudinal$sampleStatus==mysampletype,i])})
  return(mychart)
}
#adds comparisons of mean/median between two sample types, i.e. drain_1_Y and drain_1_N to see if a feature is higher in drain collection 1 in breasts that go on to be infected compared to breasts that remain uninfected
#get differences of feature mean, median (in both log10 and log2, also fold-change not in a log scale) and also P values and AUC between sample types
addcomparison <- function(mysampletype1, mysampletype2, mychart){ #difference is type1 - type2
  mychart[paste0(mysampletype1,"__",mysampletype2,"__med_diff_log10")] <- sapply(mychart$molecule, function(i){median(instudyLongitudinal[instudyLongitudinal$sampleStatus==mysampletype1,i])}) - sapply(mychart$molecule, function(i){median(instudyLongitudinal[instudyLongitudinal$sampleStatus==mysampletype2,i])})
  mychart[paste0(mysampletype1,"__",mysampletype2,"__med_diff_log2")] <- mychart[paste0(mysampletype1,"__",mysampletype2,"__med_diff_log10")]/log10(2)
  mychart[paste0(mysampletype1,"__",mysampletype2,"__mean_diff_log10")] <- sapply(mychart$molecule, function(i){mean(instudyLongitudinal[instudyLongitudinal$sampleStatus==mysampletype1,i])}) - sapply(mychart$molecule, function(i){mean(instudyLongitudinal[instudyLongitudinal$sampleStatus==mysampletype2,i])})
  mychart[paste0(mysampletype1,"__",mysampletype2,"__mean_diff_log2")] <- mychart[paste0(mysampletype1,"__",mysampletype2,"__mean_diff_log10")]/log10(2)
  mychart[paste0(mysampletype1,"__",mysampletype2,"__mean_diff_foldchange")] <- 2^mychart[paste0(mysampletype1,"__",mysampletype2,"__mean_diff_log2")]
  mychart[paste0(mysampletype1,"__",mysampletype2,"__pval_wilcoxon")] <- sapply(mychart$molecule, function(i){wilcox.test(unlist(instudyLongitudinal[instudyLongitudinal$sampleStatus==mysampletype1,i]),unlist(instudyLongitudinal[instudyLongitudinal$sampleStatus==mysampletype2,i]))$p.value})
  mychart[paste0(mysampletype1,"__",mysampletype2,"__auc")] <- sapply(mychart$molecule, function(i){roc(unlist(instudyLongitudinal[instudyLongitudinal$sampleStatus==mysampletype1 | instudyLongitudinal$sampleStatus==mysampletype2,"sampleStatus"]),
                                                                                                        unlist(instudyLongitudinal[instudyLongitudinal$sampleStatus==mysampletype1 | instudyLongitudinal$sampleStatus==mysampletype2,i]),
                                                                                                        direction="<", 
                                                                                                        levels=c(mysampletype2,mysampletype1))$auc})
  
  return(mychart)
}
#molecules of interest
longitudinalchart <- data.frame(c("AD1","pos_1685.23651_10.418","pos_215.03473_5.386","pos_461.33469_13.882","pos_143.11832_1.314","pos_244.22611_1.022","pos_333.13215_6.697","pos_218.12658_4.867","pos_216.14734_4.603","pos_556.32152_5.537","pos_1077.58075_5.559"))
colnames(longitudinalchart) <- "molecule"

#add all means/medians for future comparisons
longitudinalchart <- addmeanmed("drain_1_N",longitudinalchart)
longitudinalchart <- addmeanmed("drain_1_Y",longitudinalchart)
longitudinalchart <- addmeanmed("drain_2_N",longitudinalchart)
longitudinalchart <- addmeanmed("drain_2_Y",longitudinalchart)
longitudinalchart <- addmeanmed("drain_3_N",longitudinalchart)
longitudinalchart <- addmeanmed("drain_3_Y",longitudinalchart)
longitudinalchart <- addmeanmed("N",longitudinalchart)
longitudinalchart <- addmeanmed("Y",longitudinalchart)

#compare the progression of breasts with TEs removed due to infection vs breasts with TEs removed for non-infectious reasons at drain 1, drain 2, drain 3, and time of implant removal
longitudinalchart <- addcomparison("drain_1_Y","drain_1_N",longitudinalchart)
longitudinalchart <- addcomparison("drain_2_Y","drain_2_N",longitudinalchart)
longitudinalchart <- addcomparison("drain_3_Y","drain_3_N",longitudinalchart)
longitudinalchart <- addcomparison("Y","N",longitudinalchart)

#a quick graph showing that in general, the predictive power of drains as measured by AUCs increases over time (i.e. the later the collection and closer to infection, the better the predictive power. Which is to be expected.)
aucchart <- longitudinalchart[,c(1,grep("auc",colnames(longitudinalchart)))]
aucchart["col"] <- c("red","red4","lightgreen","darkblue","skyblue","skyblue3","orange1","orange3","orange4","cyan","pink")
plot(x=0,y=0,type="n",xlab="collection #",ylab="AUC",xlim=c(0.7,4.3),ylim=c(0.5,1))
apply(aucchart, 1, function(i){lines(c(1:4),i[2:5], col=i[6])})
legend("topleft",legend=aucchart$molecule, fill=aucchart$col, cex=0.5)

pvalchart <- longitudinalchart[,c(1,grep("pval",colnames(longitudinalchart)))]
forfdrcalc <- reshape2::melt(pvalchart, id.vars="molecule")
forfdrcalc["fdr"] <- p.adjust(forfdrcalc$value, method="BH")
forfdrcalc["fdr_colname"] <- str_replace(forfdrcalc$variable,"pval","fdr")

for(fdrname in unique(forfdrcalc$fdr_colname)){
  longitudinalchart[fdrname] <- sapply(longitudinalchart$molecule, function(molec){forfdrcalc$fdr[forfdrcalc$molecule==molec & forfdrcalc$fdr_colname==fdrname]})
}

write_csv(longitudinalchart, "figures_and_tables/table_4_longitudinalMetabolites_FULL.csv")

#make a smaller chart with only pertinent comparisons to be placed in the actual paper
longitudinalchartmini <- longitudinalchart[,c(
  "molecule",
  "drain_1_Y__drain_1_N__mean_diff_log2",
  "drain_1_Y__drain_1_N__mean_diff_foldchange",
  "drain_1_Y__drain_1_N__auc",
  "drain_1_Y__drain_1_N__pval_wilcoxon",
  "drain_1_Y__drain_1_N__fdr_wilcoxon",
  
  "drain_2_Y__drain_2_N__mean_diff_log2",
  "drain_2_Y__drain_2_N__mean_diff_foldchange",
  "drain_2_Y__drain_2_N__auc",
  "drain_2_Y__drain_2_N__pval_wilcoxon",
  "drain_2_Y__drain_2_N__fdr_wilcoxon",
  
  "drain_3_Y__drain_3_N__mean_diff_log2",
  "drain_3_Y__drain_3_N__mean_diff_foldchange",
  "drain_3_Y__drain_3_N__auc",
  "drain_3_Y__drain_3_N__pval_wilcoxon",
  "drain_3_Y__drain_3_N__fdr_wilcoxon",
  
  "Y__N__mean_diff_log2",
  "Y__N__mean_diff_foldchange",
  "Y__N__auc",
  "Y__N__pval_wilcoxon",
  "Y__N__fdr_wilcoxon"
  
)]

#rounding and/or sig figs to make things look pretty
longitudinalchartmini[,grep("log2", colnames(longitudinalchartmini))] <- round(longitudinalchartmini[,grep("log2", colnames(longitudinalchartmini))],1)
longitudinalchartmini[,grep("foldchange", colnames(longitudinalchartmini))] <- round(longitudinalchartmini[,grep("foldchange", colnames(longitudinalchartmini))],1)
longitudinalchartmini[,grep("auc", colnames(longitudinalchartmini))] <- round(longitudinalchartmini[,grep("auc", colnames(longitudinalchartmini))],2)
longitudinalchartmini[,grep("pval", colnames(longitudinalchartmini))] <- signif(longitudinalchartmini[,grep("pval", colnames(longitudinalchartmini))],2)
longitudinalchartmini[,grep("fdr", colnames(longitudinalchartmini))] <- signif(longitudinalchartmini[,grep("fdr", colnames(longitudinalchartmini))],2)
write_csv(longitudinalchartmini, "figures_and_tables/table_4_longitudinalMetabolites_SHORT.csv")
####MAKING FIGURE 6##########

#this plots feature abundance in a dot plot and boxplot, but includes drain samples as well
makeLongitudinalDotplot <- function(CID, mydf=instudyLongitudinal, mytitle=CID,fdrvec="notCalculated", myblanks=posblanks){
  comparisonsList <- list(c("drain_1_N","drain_1_Y"),c("drain_2_N","drain_2_Y"),c("drain_3_N","drain_3_Y"),c("N","Y"))#,c("drain_1_Y","drain_2_Y"),c("drain_1_Y","drain_3_Y"))
  minval <- min(mydf[,CID])
  maxval <- max(mydf[,CID])
  rangeCID <- abs(diff(range(mydf[,CID])))
  
  p2<-ggplot(mydf, aes(x=sampleStatus, y=.data[[CID]])) + 
    mytheme+
    geom_boxplot(fill=NA, color=rep(c("blue","red"),4), outliers=FALSE, linewidth=0.7)+
    geom_dotplot(binaxis='y', stackdir='center',binwidth=1/50*abs(diff(range(mydf[,CID]))))+
    ggtitle(mytitle)+ylab("log10 peak area")+theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank())+
    scale_x_discrete(labels = c("drain_1_N" = "D1\nU", "drain_1_Y" = "D1\nInf","drain_2_N" = "D2\nU","drain_2_Y" = "D2\nInf","drain_3_N" = "D3\nU","drain_3_Y" = "D3\nInf","N"="SR\nU","Y"="SR\nInf"))+
    scale_y_continuous(limits = c(minval-0.2*rangeCID, maxval+0.3*rangeCID), breaks = seq(floor(min(mydf[,CID])-0.5),ceiling(max(mydf[,CID])+0.5), by = 1))+
    geom_hline(yintercept=mean(myblanks[,CID]), linetype='dashed', col = "black")+
    geom_signif(
      comparisons = comparisonsList,
      textsize=2.8224,
      family="Helvetica",
      vjust=-0.1,
      y_position = c(max(mydf[,CID])+0.02*rangeCID, 
                     max(mydf[,CID])+0.02*rangeCID,
                     max(mydf[,CID])+0.02*rangeCID,
                     max(mydf[,CID])+0.02*rangeCID,
                     max(mydf[,CID])+0.1*rangeCID,
                     max(mydf[,CID])+0.2*rangeCID),
      map_signif_level = TRUE,
      color="black",
      annotation=fdrvec
    )
  
  plot(p2)
  return(p2)
}


getfdrsfromchart <- function(myCID){
  return(c(
    plotsignifstars(longitudinalchartmini$drain_1_Y__drain_1_N__fdr_wilcoxon[longitudinalchartmini$molecule==myCID]),
    plotsignifstars(longitudinalchartmini$drain_2_Y__drain_2_N__fdr_wilcoxon[longitudinalchartmini$molecule==myCID]),
    plotsignifstars(longitudinalchartmini$drain_3_Y__drain_3_N__fdr_wilcoxon[longitudinalchartmini$molecule==myCID]),
    plotsignifstars(longitudinalchartmini$Y__N__fdr_wilcoxon[longitudinalchartmini$molecule==myCID])
  ))
}

#"pos_218.12658_4.867","pos_216.14734_4.603"
tiff(filename = "figures_and_tables/figure_6_FINAL.tiff", width = 4, height = 5.5, units = "in", res = 600)
ggdraw() +
  draw_plot(makeLongitudinalDotplot("AD1",mytitle="HNP1", fdrvec=getfdrsfromchart("AD1")), x = 0, y = 3/4, width = 0.5, height = 1/4) +
  draw_plot(makeLongitudinalDotplot("pos_1685.23651_10.418",mytitle="HNP2", fdrvec=getfdrsfromchart("pos_1685.23651_10.418")), x = 0.5, y = 3/4, width = 0.5, height = 1/4) +
  draw_plot(makeLongitudinalDotplot("pos_215.03473_5.386",mytitle="3-chlorotyrosine", fdrvec=getfdrsfromchart("pos_215.03473_5.386")), x = 0, y = 1/2, width = 0.5, height = 1/4) +
  draw_plot(makeLongitudinalDotplot("pos_461.33469_13.882",mytitle="Hexosyl-sphinogsine", fdrvec=getfdrsfromchart("pos_461.33469_13.882")), x = 0.5, y = 1/2, width = 0.5, height = 1/4) +
  draw_plot(makeLongitudinalDotplot("pos_218.12658_4.867",mytitle="Ser-Leu", fdrvec=getfdrsfromchart("pos_218.12658_4.867")), x = 0, y = 1/4, width = 0.5, height = 1/4) +
  draw_plot(makeLongitudinalDotplot("pos_216.14734_4.603",mytitle="Val-Val", fdrvec=getfdrsfromchart("pos_216.14734_4.603")), x = 0.5, y = 1/4, width = 0.5, height = 1/4) +
  draw_plot(makeLongitudinalDotplot("pos_333.13215_6.697",mytitle="Trp-Glu", fdrvec=getfdrsfromchart("pos_333.13215_6.697")), x = 0, y = 0, width = 0.5, height = 1/4) +
  draw_plot(makeLongitudinalDotplot("pos_556.32152_5.537",mytitle="Peptide-like 556.3 Da", fdrvec=getfdrsfromchart("pos_556.32152_5.537")), x = 0.5, y = 0, width = 0.5, height = 1/4)+
  draw_plot_label(label = c("A", "B", "C","D","E","F","G","H"), size = 15,
                  x = c(0, 0.5, 0, 0.5, 0,0.5,0,0.5), y = c(1, 1, 3/4, 3/4, 1/2, 1/2, 1/4, 1/4))
dev.off()

tiff(filename = "figures_and_tables/figure_S9_FINAL.tiff", width = 2, height = 5.5*3/4, units = "in", res = 600)
ggdraw() +
  draw_plot(makeLongitudinalDotplot("pos_143.11832_1.314",mytitle="Diacetyl-spermine", fdrvec=getfdrsfromchart("pos_143.11832_1.314")), x = 0, y = 2/3, width = 1, height = 1/3) +
  draw_plot(makeLongitudinalDotplot("pos_244.22611_1.022",mytitle="Acetylspermine", fdrvec=getfdrsfromchart("pos_244.22611_1.022")), x = 0, y = 1/3, width = 1, height = 1/3) +
  draw_plot(makeLongitudinalDotplot("pos_1077.58075_5.559",mytitle="Peptide-like\n1077.6 Da", fdrvec=getfdrsfromchart("pos_1077.58075_5.559")), x = 0, y = 0, width = 1, height = 1/3) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0, 0), y = c(1, 2/3, 1/3))
dev.off()

sdv_6a <- instudyLongitudinal[,c("sampleStatus","AD1")]
sdv_6b <- instudyLongitudinal[,c("sampleStatus","pos_1685.23651_10.418")]
sdv_6c <- instudyLongitudinal[,c("sampleStatus","pos_215.03473_5.386")]
sdv_6d <- instudyLongitudinal[,c("sampleStatus","pos_461.33469_13.882")]
sdv_6e <- instudyLongitudinal[,c("sampleStatus","pos_218.12658_4.867")]
sdv_6f <- instudyLongitudinal[,c("sampleStatus","pos_216.14734_4.603")]
sdv_6g <- instudyLongitudinal[,c("sampleStatus","pos_333.13215_6.697")]
sdv_6h <- instudyLongitudinal[,c("sampleStatus","pos_556.32152_5.537")]


write.csv(sdv_6a, "supporting_data_values/sdv_6a.csv")
write.csv(sdv_6b, "supporting_data_values/sdv_6b.csv")
write.csv(sdv_6c, "supporting_data_values/sdv_6c.csv")
write.csv(sdv_6d, "supporting_data_values/sdv_6d.csv")
write.csv(sdv_6e, "supporting_data_values/sdv_6e.csv")
write.csv(sdv_6f, "supporting_data_values/sdv_6f.csv")
write.csv(sdv_6g, "supporting_data_values/sdv_6g.csv")
write.csv(sdv_6h, "supporting_data_values/sdv_6h.csv")

sdv_s9a <- instudyLongitudinal[,c("sampleStatus","pos_143.11832_1.314")]
sdv_s9b <- instudyLongitudinal[,c("sampleStatus","pos_244.22611_1.022")]
sdv_s9c <- instudyLongitudinal[,c("sampleStatus","pos_1077.58075_5.559")]

write.csv(sdv_s9a, "supporting_data_values/sdv_s9a.csv")
write.csv(sdv_s9b, "supporting_data_values/sdv_s9b.csv")
write.csv(sdv_s9c, "supporting_data_values/sdv_s9c.csv")
############################################################################################################
#FIGURE 7

#new data object so it can be altered
instudyPseud <- instudyLongitudinal
instudyPseud["xbrID"] <- unlist(rownames(instudyPseud))
#add a row on whether Pseudomonas aeruginosa was isolated from a sample
xbrkey["Pseudomonas aeruginosa"] <- sapply(xbrkey$bacteriaIsolate_ALL, function(q){grepl("Pseudomonas aeruginosa",q)})
instudyPseud["pseudoStatus"] <- sapply(instudyPseud$xbrID, function(i){paste0(xbrkey$sampletype[xbrkey$runID==i],"_",ifelse(xbrkey$`Pseudomonas aeruginosa`[xbrkey$runID==i],"PA+","PA-"))})
instudyPseud$pseudoStatus <- factor(instudyPseud$pseudoStatus, levels = c("drain_PA-","drain_PA+","seroma_PA-","seroma_PA+"))


makepseudoplot <- function(CID,myblanks=posblanks, pseudname=CID, mydf=instudyPseud){
  mydf <- instudyPseud
  myPval <- "notCalculated"
  minval <- min(mydf[,CID])
  maxval <- max(mydf[,CID])
  rangeCID <- abs(diff(range(mydf[,CID])))
  myname <- pseudname
  
  p2<-ggplot(mydf, aes(x=pseudoStatus, y=.data[[CID]])) + 
    mytheme+
    geom_boxplot(fill=NA, color=c("grey40","lightseagreen","grey40","lightseagreen"), outliers=FALSE, linewidth=0.7)+
    geom_dotplot(binaxis='y', stackdir='center',binwidth=1/80*abs(diff(range(mydf[,CID]))))+
    ggtitle(myname)+ylab("log10 peak area")+theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank())+
    scale_y_continuous(limits = c(minval-0.2*rangeCID, maxval+0.4*rangeCID), breaks = seq(floor(min(mydf[,CID])-0.5),ceiling(max(mydf[,CID])+0.5), by = 1))+
    scale_x_discrete(labels = c("drain_PA-" = "D\nPA-","drain_PA+" = "D\nPA+","seroma_PA-" = "S\nPA-", "seroma_PA+"="S\nPA+"))+
    #theme(aspect.ratio = 2/1, plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))+
    geom_hline(yintercept=mean(myblanks[,CID]), linetype='dotted', col = "black")+
    geom_signif(
      comparisons = list(c("drain_PA-", "drain_PA+"),c("seroma_PA-", "seroma_PA+")),#,c("drain_PA+","seroma_PA+")),
      map_signif_level = TRUE,
      textsize=2.8224,
      family="Helvetica",
      y_position = c(max(mydf[,CID])+0.02*rangeCID, 
                     max(mydf[,CID])+0.02*rangeCID),
      #max(mydf[,CID])+0.2*rangeCID),
      color="black",
      annotation=c(
        plotsignifstars(wilcox.test(mydf[mydf$pseudoStatus=="drain_PA-",CID],mydf[mydf$pseudoStatus=="drain_PA+",CID])$p.value),
        plotsignifstars(wilcox.test(mydf[mydf$pseudoStatus=="seroma_PA-",CID],mydf[mydf$pseudoStatus=="seroma_PA+",CID])$p.value)
        #paste0("p=",signif(wilcox.test(mydf[mydf$pseudoStatus=="drain_PA-",CID],mydf[mydf$pseudoStatus=="drain_PA+",CID])$p.value,2)),
        #paste0("p=",signif(wilcox.test(mydf[mydf$pseudoStatus=="seroma_PA-",CID],mydf[mydf$pseudoStatus=="seroma_PA+",CID])$p.value,2))#,
        #paste0("p=",signif(wilcox.test(mydf[mydf$pseudoStatus=="drain_PA+",CID],mydf[mydf$pseudoStatus=="seroma_PA+",CID])$p.value,2))
      )) # Add pairwise comparisons p-value
  
  plot(p2)
}
finalpseudomonalCIDs <- c("pos_650.38711_16.274","pos_504.32958_16.273","pos_324.06017_13.565")
finalpseudnames <- c("Di-rhamnolipid","Mono-rhamnolipid","Pyochelin")

pseudarray <- lapply(1:3, function(i){makepseudoplot(finalpseudomonalCIDs[i], pseudname = finalpseudnames[i])})

tiff("figures_and_tables/figure_7_FINAL.tiff", width=4.5,height=1.5, units="in",res=600)
ggdraw() +
  draw_plot(pseudarray[[1]], x = 0, y = 0, width = 1/3, height = 1) +
  draw_plot(pseudarray[[2]], x = 1/3, y = 0, width = 1/3, height = 1) +
  draw_plot(pseudarray[[3]], x = 2/3, y = 0, width = 1/3, height = 1) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 1/3, 2/3), y = c(1, 1, 1))
dev.off()

sdv_7a <- instudyPseud[,c("pseudoStatus","pos_650.38711_16.274")]
sdv_7b <- instudyPseud[,c("pseudoStatus","pos_504.32958_16.273")]
sdv_7c <- instudyPseud[,c("pseudoStatus","pos_324.06017_13.565")]

write.csv(sdv_7a, "supporting_data_values/sdv_7a.csv")
write.csv(sdv_7b, "supporting_data_values/sdv_7b.csv")
write.csv(sdv_7c, "supporting_data_values/sdv_7c.csv")
###################################################################################################################################
#FIGURE S5,6,8
#Cefazolin
#requires some components made above

#for the 30 negative correlates of infection at time of implant removal, see which features are correlated across samples
negcormat <- matrix(0,nrow=length(negcorCIDs),ncol=length(negcorCIDs))
colnames(negcormat) <- negcorCIDs
rownames(negcormat) <- negcorCIDs
for(i in negcorCIDs){
  for(j in negcorCIDs){
    negcormat[i,j] <- cor(finaldata[,i],finaldata[,j])
  }
}
#you can see that some negative correlates are also correlated with each other, which may be due to chemical relation (i.e. cefazolin and metabolites/impurities, crystal violet and de-methylated variants)
heatmap(negcormat, scale = "none")

# Generate the heatmap

#dotplot/boxplot of cefazolin
tiff(filename = "figures_and_tables/fig_S5a.tiff", width = 1.5, height = 1.5, units = "in", res = 600)
monochromeDotplotYN("pos_454.02951_9.9", myPval = plotsignifstars(aucdf$fdrMW[aucdf$cid=="pos_454.02951_9.9"]), mytitle = "Cefazolin")
dev.off()

sdv_s5a <- logareaall[,c("InfUninf","pos_454.02951_9.9")]
write.csv(sdv_s5a, "supporting_data_values/sdv_s5a.csv")

#dotplot/boxplot of crystal violet
tiff(filename = "figures_and_tables/fig_S6a.tiff", width = 1.5, height = 1.5, units = "in", res = 600)
monochromeDotplotYN("pos_371.23587_14.303", myPval = plotsignifstars(aucdf$fdrMW[aucdf$cid=="pos_371.23587_14.303"]), mytitle = "Crystal violet")
dev.off()

sdv_s6a <- logareaall[,c("InfUninf","pos_371.23587_14.303")]
write.csv(sdv_s6a, "supporting_data_values/sdv_s6a.csv")

#heatmap showing correlation between features, with only cefazolin labeled
#note I sometimes get a random error where the heatmap plots to RStudio instead of to the tiff. I am uncertain why, but re-running it normally works
tiff(filename = "figures_and_tables/fig_S5c.tiff", width = 5, height = 3, units = "in", res = 600)
pheatmap::pheatmap(negcormat, 
                   col = colorRampPalette(c("blue", "white", "red"))(1000),
                   breaks = seq(-1, 1, length.out = 1001),
                   scale = "none",
                   labels_col=rep("",length(colnames(negcormat))),
                   border_color = NA,
                   labels_row=sapply(rownames(negcormat), function(i){ifelse(i=="pos_454.02951_9.9","cefazolin","")}))
dev.off()

#heatmap of correlation between features, with only crystal violet labeled
tiff(filename = "figures_and_tables/fig_S6f.tiff", width = 5, height = 3, units = "in", res = 600)
pheatmap::pheatmap(negcormat, 
                   col = colorRampPalette(c("blue", "white", "red"))(1000),
                   breaks = seq(-1, 1, length.out = 1001),
                   scale = "none",
                   labels_col=rep("",length(colnames(negcormat))),
                   border_color = NA,
                   labels_row=sapply(rownames(negcormat), function(i){ifelse(i=="pos_371.23587_14.303","crystal violet","")}))
dev.off()

#for the group of things correlated with cefazolin, make a table with the feature IDs, mz, RT, etc
cefcorrelates <- rownames(negcormat)[negcormat["pos_454.02951_9.9",]>0.7]

rawposcompounds <- read_csv("input_files/20240116_all_implant_samples_pos_MS1_aln_rerun3_compounds.csv")
rawposcompounds["CID"] <- paste0("pos_",rawposcompounds$`Calc. MW`,"_",rawposcompounds$`RT [min]`)

rawneg <- read_csv("input_files/20240131_all_implant_samples_neg_MS2_aln_noRTalign_compounds.csv")
rawneg["CID"] <- paste0("neg_",rawneg$`Calc. MW`,"_",rawneg$`RT [min]`)

cefmat <- data.frame(cbind(cefcorrelates,
                           sapply(cefcorrelates, function(i){
                             ifelse(substr(i,1,3)=="pos",
                                    rawposcompounds$`m/z`[rawposcompounds$CID==i],
                                    rawneg$`m/z`[rawneg$CID==i])
                           }),
                           sapply(cefcorrelates, function(i){
                             ifelse(substr(i,1,3)=="pos",
                                    rawposcompounds$`Calc. MW`[rawposcompounds$CID==i],
                                    rawneg$`Calc. MW`[rawneg$CID==i])
                           }),
                           sapply(cefcorrelates, function(i){
                             ifelse(substr(i,1,3)=="pos",
                                    rawposcompounds$`RT [min]`[rawposcompounds$CID==i],
                                    rawneg$`RT [min]`[rawneg$CID==i])
                           }),
                           round(sapply(cefcorrelates, function(i){negcormat["pos_454.02951_9.9",i]}),2)
))

colnames(cefmat) <- c("feature","m/z","Calc. MW","RT","correlation to cefazolin")
cefmat <- cefmat[order(cefmat$`correlation to cefazolin`, decreasing=TRUE),]

write_csv(cefmat, "figures_and_tables/table_s1_cefazolin.csv")

#do the same for crystal violet
cvcorrelates <- rownames(negcormat)[negcormat["pos_371.23587_14.303",]>0.7]
cvmat <- data.frame(cbind(cvcorrelates,
                          sapply(cvcorrelates, function(i){
                            ifelse(substr(i,1,3)=="pos",
                                   rawposcompounds$`m/z`[rawposcompounds$CID==i],
                                   rawneg$`m/z`[rawneg$CID==i])
                          }),
                          sapply(cvcorrelates, function(i){
                            ifelse(substr(i,1,3)=="pos",
                                   rawposcompounds$`Calc. MW`[rawposcompounds$CID==i],
                                   rawneg$`Calc. MW`[rawneg$CID==i])
                          }),
                          sapply(cvcorrelates, function(i){
                            ifelse(substr(i,1,3)=="pos",
                                   rawposcompounds$`RT [min]`[rawposcompounds$CID==i],
                                   rawneg$`RT [min]`[rawneg$CID==i])
                          }),
                          round(sapply(cvcorrelates, function(i){negcormat["pos_371.23587_14.303",i]}),2)
))

colnames(cvmat) <- c("feature","m/z","Calc. MW","RT","correlation to crystal violet")
cvmat <- cvmat[order(cvmat$`correlation to crystal violet`,decreasing=TRUE),]
write_csv(cvmat, "figures_and_tables/table_s2_cv.csv")


#make a plot of bromotyrosine for the supplement
tiff(filename = "figures_and_tables/fig_S8c.tiff", width = 1.5, height = 1.5, units = "in", res = 600)
brTyrCID <- "pos_258.98419_5.868"
monochromeDotplotYN(brTyrCID,mydf=logareaallADs, mytitle="3-bromotyrosine", myPval = plotsignifstars(wilcox.test(unlist(logareaallADs[logareaallADs$InfUninf=="Infected",brTyrCID]),unlist(logareaallADs[logareaallADs$InfUninf=="Uninfected",brTyrCID]))$p.value))

dev.off()
sdv_s8c <- logareaall[,c("InfUninf","pos_258.98419_5.868")]
write.csv(sdv_s8c, "supporting_data_values/sdv_s8c.csv")
