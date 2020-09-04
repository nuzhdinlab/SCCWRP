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

#Read in morphological stream community data and format them as presence/absence tables.
AlgalInput <- read.table("AlgTaxa.csv", header=T, sep=",",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")

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

#Read in morphological stream community data and format them as presence/absence tables.
BMIInput <- read.table("BugTaxa.csv", header=T, sep=",",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")

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
