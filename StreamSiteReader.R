library("plyr")
library(dplyr)
library("ggplot2")
library(lubridate)
library("ape")
library("vegan")
library("microbiome")
library(data.table)
library(tidyr)

setwd("~/Desktop/SCCWRP")
#Read in site data containing biological counts, water chemistry, and land usage
#values.  Generate a merged data set.

#GISBiochemData <- read.table("GISBiochemData.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)

#If you need to aggregate site data please proceed here.
#Read in algae data from SMC sites.
algaeDataSMCRaw <- read.table("AlgaeTax_dnaSites_SMC.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
algaeDataSMC <- filter(algaeDataSMCRaw, Replicate==1)
#Subset columns of interest for the SMC sites.
algaeDataSMC <- algaeDataSMCRaw[,c(1,3,43,40)]
#Change the header name for station ID.
names(algaeDataSMC)[names(algaeDataSMC)=="Sample Station ID"]<-"SampleStationID"
#Determine the algal totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ SampleStationID,algaeDataSMC))
colnames(tmp) <- c("SampleStationID","ActualOrganismCount")
#Add algal totals count column to algae dataframe.
algaeDataSMC <- merge(algaeDataSMC,tmp,"SampleStationID")
#Calculate the relative abundance of algae data.
algaeDataSMC$Measurement <- with(algaeDataSMC,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
algaeDataSMC$MeasurementType <- with(algaeDataSMC,"Algal relative abundance")
#Force a uniform date format
algaeDataSMC$SampleDate <- mdy(algaeDataSMC$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
algaeDataSMC$UniqueID <- with(algaeDataSMC,paste(algaeDataSMC$SampleStationID,"SMC",algaeDataSMC$SampleDate))
#Find sampling year.
algaeDataSMC$Year <- year(algaeDataSMC$SampleDate)

#Read in algae data from CEDEN sites.
algaeDataCEDENRaw <- read.table("AlgaeTax_dnaSites_CEDEN.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
algaeDataCEDEN <- filter(algaeDataCEDENRaw, CollectionReplicate==1)
#Subset columns of interest for the CEDEN sites.
algaeDataCEDEN <- algaeDataCEDENRaw[,c(6,11,36,26)]
#Change names to uniforma schema.
names(algaeDataCEDEN)[names(algaeDataCEDEN)=="StationCode"]<-"SampleStationID"
#names(algaeDataCEDEN)[names(algaeDataCEDEN)=="Species"]<-"FinalID"
#Determine the algal totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ SampleStationID,algaeDataCEDEN))
colnames(tmp) <- c("SampleStationID","ActualOrganismCount")
#Add algal totals count column to algae dataframe.
algaeDataCEDEN <- merge(algaeDataCEDEN,tmp,"SampleStationID")
#Calculate the relative abundance of algae data.
algaeDataCEDEN$Measurement <- with(algaeDataCEDEN,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
algaeDataCEDEN$MeasurementType <- with(algaeDataCEDEN,"Algal relative abundance")
#Force a uniform date format
algaeDataCEDEN$SampleDate <- ymd(algaeDataCEDEN$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
algaeDataCEDEN$UniqueID <- with(algaeDataCEDEN,paste(algaeDataCEDEN$SampleStationID,"CEDEN",algaeDataCEDEN$SampleDate))
#Find sampling year.
algaeDataCEDEN$Year <- year(algaeDataCEDEN$SampleDate)

#The SWAMP data file is in a somewhat irregular format and this is accounted for
#when being read in.
algaeDataSWAMPRaw <- read.table("AlgaeTaxonomy_dnaSamples_SWAMP.csv", fill=TRUE,header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
algaeDataSWAMP <- filter(algaeDataSWAMPRaw, Replicate==1)
algaeDataSWAMP <- algaeDataSWAMPRaw[,c(6,8,97,90)]
#Change names to uniforma schema.
names(algaeDataSWAMP)[names(algaeDataSWAMP)=="StationCode"]<-"SampleStationID"
#Determine the algal totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ SampleStationID,algaeDataSWAMP))
colnames(tmp) <- c("SampleStationID","ActualOrganismCount")
#Add algal totals count column to algae dataframe.
algaeDataSWAMP <- merge(algaeDataSWAMP,tmp,"SampleStationID")
#Calculate the relative abundance of algae data.
algaeDataSWAMP$Measurement <- with(algaeDataSWAMP,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
algaeDataSWAMP$MeasurementType <- with(algaeDataSWAMP,"Algal relative abundance")
#Force a uniform date format
algaeDataSWAMP$SampleDate <- mdy(algaeDataSWAMP$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
algaeDataSWAMP$UniqueID <- with(algaeDataSWAMP,paste(algaeDataSWAMP$SampleStationID,"SWAMP",algaeDataSWAMP$SampleDate))
#Find sampling year.
algaeDataSWAMP$Year <- year(algaeDataSWAMP$SampleDate)

#Create merged algae data set.
algaeData <- do.call("rbind",list(algaeDataSMC,algaeDataSWAMP,algaeDataCEDEN))
algaeData <- na.omit(algaeData)
#Determine the insect Shannon diversity and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(-Measurement*log(Measurement) ~ UniqueID,algaeData))
colnames(tmp) <- c("UniqueID","Shannon")
algaeData <- merge(algaeData,tmp,"UniqueID")

#Read in insect data from SMC sites.
insectDataSMCRAW <- read.csv("BugTax_dnaSites_SMC.csv")
#Subset only replicate 1
insectDataSMC <- filter(insectDataSMCRAW, FieldReplicate==1)
#Subset columns of interest.
insectDataSMC <- insectDataSMCRAW[,c(1,3,9,6)]
#Change names to uniforma schema.
names(insectDataSMC)[names(insectDataSMC)=="Sample.Station.ID"]<-"SampleStationID"
#Determine the insect totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ SampleStationID,insectDataSMC))
colnames(tmp) <- c("SampleStationID","ActualOrganismCount")
#Add insect totals count column to insect dataframe.
insectDataSMC <- merge(insectDataSMC,tmp,"SampleStationID")
#Calculate the relative abundance of insect data.
insectDataSMC$Measurement <- with(insectDataSMC,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
insectDataSMC$MeasurementType <- with(insectDataSMC,"Invertebrate relative abundances")
#Force a uniform date format
insectDataSMC$SampleDate <- mdy(insectDataSMC$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
insectDataSMC$UniqueID <- with(insectDataSMC,paste(insectDataSMC$SampleStationID,"SMC",insectDataSMC$SampleDate))
#Find sampling year.
insectDataSMC$Year <- year(insectDataSMC$SampleDate)

#Read in insect data from CEDEN sites.
insectDataCEDENRAW <- read.table("BugTax_dnaSites_CEDEN.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
insectDataCEDEN <- filter(insectDataCEDENRAW, CollectionReplicate==1)
#Subset columns of interest.
insectDataCEDEN <- insectDataCEDENRAW[,c(6,11,36,26)]
#Change names to uniforma schema.
names(insectDataCEDEN)[names(insectDataCEDEN)=="StationCode"]<-"SampleStationID"
#Determine the insect totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ SampleStationID,insectDataCEDEN))
colnames(tmp) <- c("SampleStationID","ActualOrganismCount")
#Add insect totals count column to insect dataframe.
insectDataCEDEN <- merge(insectDataCEDEN,tmp,"SampleStationID")
#Calculate the relative abundance of insect data.
insectDataCEDEN$Measurement <- with(insectDataCEDEN,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
insectDataCEDEN$MeasurementType <- with(insectDataCEDEN,"Invertebrate relative abundance")
#Force a uniform date format
insectDataCEDEN$SampleDate <- ymd(insectDataCEDEN$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
insectDataCEDEN$UniqueID <- with(insectDataCEDEN,paste(insectDataCEDEN$SampleStationID,"CEDEN",insectDataCEDEN$SampleDate))
#Find sampling year.
insectDataCEDEN$Year <- year(insectDataCEDEN$SampleDate)

#Read in insect data from SWAMP sites.
insectDataSWAMPRAW <- read.table("BugTaxonomy_dnaSamples_SWAMP.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Subset only replicate 1
insectDataSWAMP <- filter(insectDataSWAMPRAW, Replicate==1)
#Subset columns of interest.
insectDataSWAMP <- insectDataSWAMPRAW[,c(1,8,97,90)]
#Change names to uniforma schema.
names(insectDataSWAMP)[names(insectDataSWAMP)=="Sample Station ID"]<-"SampleStationID"
#Determine the insect totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ SampleStationID,insectDataSWAMP))
colnames(tmp) <- c("SampleStationID","ActualOrganismCount")
#Add insect totals count column to insect dataframe.
insectDataSWAMP <- merge(insectDataSWAMP,tmp,"SampleStationID")
#Calculate the relative abundance of insect data.
insectDataSWAMP$Measurement <- with(insectDataSWAMP,BAResult/ActualOrganismCount)
#Add organism type for later use in merged data sets.
insectDataSWAMP$MeasurementType <- with(insectDataSWAMP,"Invertebrate relative abundance")
#Force a uniform date format
insectDataSWAMP$SampleDate <- mdy(insectDataSWAMP$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
insectDataSWAMP$UniqueID <- with(insectDataSWAMP,paste(insectDataSWAMP$SampleStationID,"SWAMP",insectDataSWAMP$SampleDate))
#Find sampling year.
insectDataSWAMP$Year <- year(insectDataSWAMP$SampleDate)

#Create merged insect data set.
insectData <- do.call("rbind",list(insectDataSMC,insectDataSWAMP,insectDataCEDEN))
insectData <- na.omit(insectData)
insectData$UniqueTaxa <- paste(insectData$UniqueID,insectData$FinalID)
#Determine the insect totals count column and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(BAResult ~ UniqueTaxa,insectData))
colnames(tmp) <- c("UniqueTaxa","BAResult")
#Add insect totals count column to insect dataframe.
insectData <- merge(subset(insectData,select=-c(3)),tmp,"UniqueTaxa")
#Calculate the relative abundance of insect data.
insectData$Measurement <- with(insectData,BAResult/ActualOrganismCount)
insectData <- insectData[!duplicated(insectData),]
#Remove UniqueTaxa column
insectData <- subset(insectData,select=-c(1))
#Reorder columns prior to merger.
insectData <- insectData[c("SampleStationID","SampleDate","BAResult","FinalID","ActualOrganismCount","Measurement","MeasurementType","UniqueID","Year")]
#Determine the insect Shannon diversity and make it a temporary dataframe.
tmp <- as.data.frame(xtabs(-Measurement*log(Measurement) ~ UniqueID,insectData))
colnames(tmp) <- c("UniqueID","Shannon")
insectData <- merge(insectData,tmp,"UniqueID")

#Merge insect and algae data.
bioData <- do.call("rbind",list(insectData,algaeData))
#bioData <- bioData[,-c(3,5)]
bioData <- bioData[!duplicated(bioData),]

#Read in chemical data for the SMC test sites.
chemDataSMCRAW <- read.table("Chem_dnaSites_SMC.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset only replicate 1
chemDataSMC <- filter(chemDataSMCRAW, FieldReplicate==1)
#Subset columns.
chemDataSMC <- chemDataSMCRAW[,c(1,3,13,17,15)]
#Introduce common naming schema.
names(chemDataSMC)[names(chemDataSMC)=="Sample Station ID"]<-"SampleStationID"
names(chemDataSMC)[names(chemDataSMC)=="AnalyteName"]<-"FinalID"
names(chemDataSMC)[names(chemDataSMC)=="Result"]<-"Measurement"
names(chemDataSMC)[names(chemDataSMC)=="Unit"]<-"MeasurementType"
#Force a uniform date format
chemDataSMC$SampleDate <- mdy(chemDataSMC$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
chemDataSMC$UniqueID <- with(chemDataSMC,paste(chemDataSMC$SampleStationID,"SMC",chemDataSMC$SampleDate))
#Find sampling year.
chemDataSMC$Year <- year(chemDataSMC$SampleDate)

#Read in chemical data for the CEDEN test sites.
chemDataCEDENRAW <- read.csv("Chem_dnaSites_CEDEN.csv")
#Subset only replicate 1
chemDataCEDEN <- filter(chemDataCEDENRAW, CollectionReplicate==1)
#Subset the data.
chemDataCEDEN <- chemDataCEDEN[,c(5,6,18,20,19)]
#Introduce common naming schema.
names(chemDataCEDEN)[names(chemDataCEDEN)=="StationCode"]<-"SampleStationID"
names(chemDataCEDEN)[names(chemDataCEDEN)=="Analyte"]<-"FinalID"
names(chemDataCEDEN)[names(chemDataCEDEN)=="Result"]<-"Measurement"
names(chemDataCEDEN)[names(chemDataCEDEN)=="Unit"]<-"MeasurementType"
#Force a uniform date format
chemDataCEDEN$SampleDate <- as.Date(chemDataCEDEN$SampleDate,format="%m/%d/%y")
#Create unique ID combining the sample station ID and sampling date.
chemDataCEDEN$UniqueID <- with(chemDataCEDEN,paste(chemDataCEDEN$SampleStationID,"CEDEN",chemDataCEDEN$SampleDate))
#Find sampling year.
chemDataCEDEN$Year <- year(chemDataCEDEN$SampleDate)

#Read in chemical data for the SWAMP test sites.
chemDataSWAMPRAW <- read.csv("Chem_dnaSamples_SWAMP.csv")
#Subset only replicate 1
chemDataSWAMP <- filter(chemDataSWAMPRAW,Replicate==1)
#Subset the data.
chemDataSWAMP <- chemDataSWAMP[,c(1,8,70,89,74)]
#Introduce common naming schema.
names(chemDataSWAMP)[names(chemDataSWAMP)=="Sample.Station.ID"]<-"SampleStationID"
names(chemDataSWAMP)[names(chemDataSWAMP)=="AnalyteName"]<-"FinalID"
names(chemDataSWAMP)[names(chemDataSWAMP)=="Result"]<-"Measurement"
names(chemDataSWAMP)[names(chemDataSWAMP)=="UnitName"]<-"MeasurementType"
#Force a uniform date format
chemDataSWAMP$SampleDate <- mdy(chemDataSWAMP$SampleDate)
#Create unique ID combining the sample station ID and sampling date.
chemDataSWAMP$UniqueID <- with(chemDataSWAMP,paste(chemDataSWAMP$SampleStationID,"SWAMP",chemDataSWAMP$SampleDate))
#Find sampling year.
chemDataSWAMP$Year <- year(chemDataSWAMP$SampleDate)

#Create merged chemical data frame.
chemData <- do.call("rbind",list(chemDataSMC,chemDataSWAMP,chemDataCEDEN))
chemData <- na.omit(chemData)

#Merge chemical and biological data.
bioChemData <- do.call("rbind",list(bioData,chemData))
#Sort bio-chem data by year.
bioChemData <- bioChemData[order(as.numeric(bioChemData$Year)),]

#Read in geospatial data.
GISDataRAW <- read.table("GIS_dnaSites.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
#Subset columns of interest.
GISData <- GISDataRAW[,-c(1,3:5,8:10,15)]
names(GISData)[names(GISData)=="StationCode"]<-"SampleStationID"
names(GISData)[names(GISData)=="New_Lat"]<-"Latitude"
names(GISData)[names(GISData)=="New_Long"]<-"Longitude"

#Merge geospatial data with biological observations.
GISBiochemData <- join(bioChemData,GISData,by="SampleStationID")
#GISBiochemData <- GISBiochemData[,-c(10:11,14:22,47:51,82:90,100)]
#Sort merged data set by year then measurement name.
GISBiochemData <- as.data.frame(GISBiochemData[order(as.numeric(GISBiochemData$Year),as.character(GISBiochemData$FinalID)),])

#Filter out erroneous negative values for physical parameter data.
chemID <- unique(chemData$FinalID)
GISBiochemData <- filter(GISBiochemData,Measurement>=0)

#Filter out duplicate rows.
GISBiochemData <- GISBiochemData[!duplicated(GISBiochemData),]

#Incorporate CSCI data.
csciData1 <- read.csv("CSCI_dnaSites.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
csciData1 <- filter(csciData1,CSCI!="NA")
csciData1 <- subset(csciData1,select=-c(StationCode))
#Fix date format.
csciData1$SampleDate <- mdy(csciData1$SampleDate)
#Create Year column.
csciData1$Year <- year(csciData1$SampleDate)
#Insert Unique ID
csciData1$UniqueID <- with(csciData1,paste(csciData1$SampleStationID,"SWAMP",csciData1$SampleDate))
#Subset UniqueID and CSCI.
csciData1 <- csciData1[,c("UniqueID","CSCI")]
#Incorporate CSCI data.
csciData2 <- read.csv("csci_scored_sites_tbl.csv", header=TRUE, sep=",",as.is=T,check.names=FALSE)
csciData2 <- filter(csciData2,CSCI!="NA")
names(csciData2)[names(csciData2)=="StationCode"]<-"SampleStationID"
#Fix date format.
csciData2$SampleDate <- mdy(csciData2$SAMPLEDATE)
#Create Year column.
csciData2$Year <- year(csciData2$SampleDate)
#Insert Unique ID
csciData2$UniqueID <- with(csciData2,paste(csciData2$SampleStationID,"SWAMP",csciData2$SampleDate))
#Subset UniqueID and CSCI.
csciData2 <- csciData2[,c("UniqueID","CSCI")]
#Merge into a single CSCI dataframe.
csciData <- rbind(csciData1,csciData2)
csciData <- csciData[!duplicated(csciData),]

#Merge CSCI data by UniqueID
GISBiochemData <- join(GISBiochemData,csciData,by=c("UniqueID"))
GISBiochemData <- GISBiochemData[!duplicated(GISBiochemData),]

#Calculate land usage index based on 1K, 5K, and catchment zone values.
#Use land usage data from 2011.
GISBiochemData$LU_2011_1K <- with(GISBiochemData,Ag_2011_1K+CODE_21_2011_1K+URBAN_2011_1K)
GISBiochemData$LU_2011_5K <- with(GISBiochemData,Ag_2011_5K+CODE_21_2011_5K+URBAN_2011_5K)
GISBiochemData$LU_2011_WS <- with(GISBiochemData,Ag_2011_WS+CODE_21_2011_WS+URBAN_2011_WS)

#Filter out data without a land usage index.
GISBiochemData <- subset(GISBiochemData,LU_2011_1K!="NA")
GISBiochemData <- subset(GISBiochemData,LU_2011_5K!="NA")
GISBiochemData <- subset(GISBiochemData,LU_2011_WS!="NA")

#Subset site data based on land usage index within 1K catchment zones.
#LD = low disturbance.  Land usage index is less than 5%.
GISBiochemDataLD1K <- GISBiochemData[which(GISBiochemData$LU_2011_1K < 5),]
GISBiochemDataLD1K$LUCategory <- "LD1K"
#MD = low disturbance.  Land usage index is between 5% and 15%.
GISBiochemDataMD1K <- GISBiochemData[which(GISBiochemData$LU_2011_1K < 15 & GISBiochemData$LU_2011_1K >= 5),]
GISBiochemDataMD1K$LUCategory <- "MD1K"
#HD = low disturbance.  Land usage index is greater than 15%.
GISBiochemDataHD1K <- GISBiochemData[which(GISBiochemData$LU_2011_1K >= 15),]
GISBiochemDataHD1K$LUCategory <- "HD1K"

#Subset site data based on land usage index within 5K catchment zones.
#LD = low disturbance.  Land usage index is less than 5%.
GISBiochemDataLD5K <- GISBiochemData[which(GISBiochemData$LU_2011_5K < 5),]
GISBiochemDataLD5K$LUCategory <- "LD5K"
#MD = low disturbance.  Land usage index is between 5% and 15%.
GISBiochemDataMD5K <- GISBiochemData[which(GISBiochemData$LU_2011_5K < 15 & GISBiochemData$LU_2011_5K >= 5),]
GISBiochemDataMD5K$LUCategory <- "MD5K"
#HD = low disturbance.  Land usage index is greater than 15%.
GISBiochemDataHD5K <- GISBiochemData[which(GISBiochemData$LU_2011_5K >= 15),]
GISBiochemDataHD5K$LUCategory <- "HD5K"

#Subset site data based on land usage index within the full water drainage catchment zones.
#LD = low disturbance.  Land usage index is less than 5%.
GISBiochemDataLDWS <- GISBiochemData[which(GISBiochemData$LU_2011_WS < 5),]
GISBiochemDataLDWS$LUCategory <- "LDWS"
#MD = low disturbance.  Land usage index is between 5% and 15%.
GISBiochemDataMDWS <- GISBiochemData[which(GISBiochemData$LU_2011_WS < 15 & GISBiochemData$LU_2011_WS >= 5),]
GISBiochemDataMDWS$LUCategory <- "MDWS"
#HD = low disturbance.  Land usage index is greater than 15%.
GISBiochemDataHDWS <- GISBiochemData[which(GISBiochemData$LU_2011_WS >= 15),]
GISBiochemDataHDWS$LUCategory <- "HDWS"

#Merge land usage subsets back together for later analytical tools.
GISBiochemData <- do.call("rbind",list(GISBiochemDataLD1K,GISBiochemDataMD1K,GISBiochemDataHD1K,GISBiochemDataLD5K,GISBiochemDataMD5K,GISBiochemDataHD5K,GISBiochemDataLDWS,GISBiochemDataMDWS,GISBiochemDataHDWS))

#Write out merged data set to read back in the future as opposed to 
#generating it each time from component data files.
write.csv(GISBiochemData,file="GISBiochemData.csv",row.names=FALSE)
