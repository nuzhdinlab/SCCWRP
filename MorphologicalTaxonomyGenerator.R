rm(list=ls())
require("plyr")
require(dplyr)
require(stringr)
require(tidyr)
require(naniar)
require(taxize)

options(ENTREZ_KEY="3f9e29c7c03842f72cf7523e34390d9f2208")
wd <- "~/Desktop/SCCWRP/MorphologicalBMIsAlgae/"
setwd(wd)
#This script generates full taxonomies for morphologically classified stream communities from SCCWRP.

#Read in morphological stream community data describing all algae (soft algae and diatoms)
#and format them as presence/absence tables.
AlgalInput <- read.table("AlgTaxa.csv", header=T, sep=",",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")
#Standarize date format.
AlgalInput$sampledate <- as.Date(AlgalInput$sampledate,format="%Y-%m-%d")
AlgalInput$sampledate <- format(AlgalInput$sampledate,format="%m/%d/%y")
#Create unique sample identifier.
AlgalInput$UniqueID <- paste(AlgalInput$stationcode,AlgalInput$sampledate)

#Read in table linking sample IDs in the metagenomic table to sample station codes.
#This will filter the morphological data to containing the same sample sets.
sampleIDs <- read.table("SampleStationCodesID.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
#Remove dubious samples
sampleIDs <- sampleIDs[which(sampleIDs$SampleNum!=202 & sampleIDs$SampleNum!=230),]
sampleIDs$Date <- as.Date(sampleIDs$Date,format="%m/%d/%y")
sampleIDs$Date <- format(sampleIDs$Date,format="%m/%d/%y")
#Create unique sample identifier.
sampleIDs$UniqueID <- paste(sampleIDs$StationCode,sampleIDs$Date)

#Filter morphological algal data to containing the same samples as the metagenomic ones.
AlgalInput <- AlgalInput[AlgalInput$UniqueID %in% sampleIDs$UniqueID,]
  
#Get a list of unique OTU names from count tables and convert to a data frame.
uniqueOTUs <- as.data.frame(unique(AlgalInput$finalid))
colnames(uniqueOTUs) <- c("FinalID")
uniqueOTUs <- as.data.frame(uniqueOTUs[uniqueOTUs$FinalID!="",])
colnames(uniqueOTUs) <- c("FinalID")
uniqueOTUs$FinalID <- as.character(uniqueOTUs$FinalID)

#Help prevent server timeouts.
httr::set_config(httr::config(http_version = 0))

# Get full taxonomies for all algal OTUs.
AlgalTaxonomies <- data.frame()
i=0
for(name in unique(uniqueOTUs$FinalID)){
  tmp <- classification(name,db="ncbi",rows=1)
  Sys.sleep(0.2)
  if(nrow(as.data.frame(tmp[1]))>1){
    tmp <- as.data.frame(tmp[1])
    colnames(tmp) <- c("taxa","rank","id")
    tmp <- tmp[tmp$rank!="clade",]
    tmp <- tmp[tmp$rank!="no rank",]
    if(nrow(tmp)>1){
      tmp <- as.data.frame(t(tmp[,c("taxa","rank")]))
      colnames(tmp) <- as.character(unlist(tmp["rank",]))
      tmp <- tmp[!row.names(tmp) %in% c("rank"),]
      rownames(tmp) <- name
      tmp$LeafTaxa <- name
      AlgalTaxonomies <- dplyr::bind_rows(AlgalTaxonomies,tmp)
      i=i+1
      print(paste("Taxon",i,length(unique(uniqueOTUs$FinalID))))
    } else{
      tmp <- classification(name,db="gbif",rows=1)
      Sys.sleep(0.2)
      if(nrow(as.data.frame(tmp[1]))>1){
        tmp <- as.data.frame(tmp[1])
        colnames(tmp) <- c("taxa","rank","id")
        tmp <- tmp[tmp$rank!="clade",]
        tmp <- tmp[tmp$rank!="no rank",]
        if(nrow(tmp)>1){
          tmp <- as.data.frame(t(tmp[,c("taxa","rank")]))
          colnames(tmp) <- as.character(unlist(tmp["rank",]))
          tmp <- tmp[!row.names(tmp) %in% c("rank"),]
          rownames(tmp) <- name
          tmp$LeafTaxa <- name
          AlgalTaxonomies <- dplyr::bind_rows(AlgalTaxonomies,tmp)
          i=i+1
          print(paste("Taxon",i,length(unique(uniqueOTUs$FinalID))))
        }
      }
    }
  }
}
#Format taxonomy data frame.
AlgalTaxonomies[] <- lapply(AlgalTaxonomies,as.character)

#Write out algal taxonomy data.
write.table(AlgalTaxonomies,"AlgalTaxonomiesMorphological.txt",quote=FALSE,sep="\t",row.names = FALSE)

#Generate full taxonomies for BMIs.
#Read in morphological stream community data for BMIS and format them as presence/absence tables.
BMIInput <- read.table("BugTaxa.csv", header=T, sep=",",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")

#Standarize date format.
BMIInput$sampledate <- gsub(" 00:00:00","",BMIInput$sampledate)
BMIInput$sampledate <- as.Date(BMIInput$sampledate,format="%Y-%m-%d")
BMIInput$sampledate <- format(BMIInput$sampledate,format="%m/%d/%y")
#Create unique sample identifier.
BMIInput$UniqueID <- paste(BMIInput$stationcode,BMIInput$sampledate)

#Filter morphological algal data to containing the same samples as the metagenomic ones.
BMIInput <- BMIInput[BMIInput$UniqueID %in% sampleIDs$UniqueID,]

#Get a list of unique OTU names from count tables and convert to a data frame.
uniqueOTUs <- as.data.frame(unique(BMIInput$finalid))
colnames(uniqueOTUs) <- c("FinalID")
uniqueOTUs <- as.data.frame(uniqueOTUs[uniqueOTUs$FinalID!="",])
colnames(uniqueOTUs) <- c("FinalID")
uniqueOTUs$FinalID <- as.character(uniqueOTUs$FinalID)

# Get full taxonomies for all BMI OTUs.
BMITaxonomies <- data.frame()
i=0
for(name in unique(uniqueOTUs$FinalID)){
  tmp <- classification(name,db="ncbi",rows=1)
  Sys.sleep(0.2)
  if(nrow(as.data.frame(tmp[1]))>1){
    tmp <- as.data.frame(tmp[1])
    colnames(tmp) <- c("taxa","rank","id")
    tmp <- tmp[tmp$rank!="clade",]
    tmp <- tmp[tmp$rank!="no rank",]
    if(nrow(tmp)>1){
      tmp <- as.data.frame(t(tmp[,c("taxa","rank")]))
      colnames(tmp) <- as.character(unlist(tmp["rank",]))
      tmp <- tmp[!row.names(tmp) %in% c("rank"),]
      rownames(tmp) <- name
      tmp$LeafTaxa <- name
      BMITaxonomies <- dplyr::bind_rows(BMITaxonomies,tmp)
      i=i+1
      print(paste("Taxon",i,length(unique(uniqueOTUs$FinalID))))
    } else{
      tmp <- classification(name,db="gbif",rows=1)
      Sys.sleep(0.2)
      if(nrow(as.data.frame(tmp[1]))>1){
        tmp <- as.data.frame(tmp[1])
        colnames(tmp) <- c("taxa","rank","id")
        tmp <- tmp[tmp$rank!="clade",]
        tmp <- tmp[tmp$rank!="no rank",]
        if(nrow(tmp)>1){
          tmp <- as.data.frame(t(tmp[,c("taxa","rank")]))
          colnames(tmp) <- as.character(unlist(tmp["rank",]))
          tmp <- tmp[!row.names(tmp) %in% c("rank"),]
          rownames(tmp) <- name
          tmp$LeafTaxa <- name
          BMITaxonomies <- dplyr::bind_rows(BMITaxonomies,tmp)
          i=i+1
          print(paste("Taxon",i,length(unique(uniqueOTUs$FinalID))))
        }
      }
    }
  }
}
#Format taxonomy data frame.
BMITaxonomies[] <- lapply(BMITaxonomies,as.character)

#Write out BMI taxonomy data.
write.table(BMITaxonomies,"BMITaxonomiesMorphological.txt",quote=FALSE,sep="\t",row.names = FALSE)
