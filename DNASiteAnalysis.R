require(dplyr)
require(plyr)

##Part 1: Run this to generate the input file for use in QGIS.
#This part is get a file for use in merging in HUC 8 watershed data through QGIS.
setwd("~/Desktop/SCCWRP/HUC8Watersheds")
#Read in updated SCCWRP site locations.
DNASiteLocations <- read.table("UpdatedDNASiteLocations.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Removed undefined entries
DNASiteLocations <- na.omit(DNASiteLocations)
#Create unique sample identifier.
DNASiteLocations$UniqueID <- paste(DNASiteLocations$StationCode,DNASiteLocations$Date.Collected)
DNASiteLocations <- DNASiteLocations[,c("Latitude","Longitude","UniqueID")]
#Read in site metadata.
DNASiteMetadata <- read.table("RefGISdata_algaeDNAsites.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Create unique sample identifier.
DNASiteMetadata$UniqueID <- paste(DNASiteMetadata$StationCode,DNASiteMetadata$Date.Collected)
#Remove outdated latitude and longitude
DNASiteMetadata <- DNASiteMetadata[, !names(DNASiteMetadata) %in% c("Latitude","Longitude")]
#Create merged data frame to output to QGIS.
MergedDNASiteData <- left_join(DNASiteLocations,DNASiteMetadata,by=c("UniqueID"))
#Output file to QGIS.
write.table(MergedDNASiteData,"UpdatedGISDNAsites.csv",quote=FALSE,sep=",",row.names = FALSE)

##Part 2: Merge in HUC 8 watershed data in QGIS.  Output the merged file.

##Part 3: Analyze the merged file generated via Step 1 and Step 2.
setwd("~/Desktop/SCCWRP")
#setwd("/home/cmb-07/sn1/alsimons/SCCWRP")

#Read in sample metadata.
SCCWRP <- read.table("SCCWRPDNASites.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Remove duplicate rows
SCCWRP <- SCCWRP[!duplicated(SCCWRP),]
#Create HUC 6 watershed column.
SCCWRP$HUC6 <- substr(SCCWRP$HUC8,1,nchar(SCCWRP$HUC8)-2)
#Create HUC 4 watershed column.
SCCWRP$HUC4 <- substr(SCCWRP$HUC8,1,nchar(SCCWRP$HUC8)-4)
#Create aggregated land use column.  Use 2011 land use data within 5 km clipped buffer.
SCCWRP$LU <- SCCWRP$code_21_2011_5k+SCCWRP$ag_2011_5k+SCCWRP$urban_2011_5k
#Keep entries with valid land use data.
SCCWRP <- SCCWRP[!is.na(SCCWRP$LU),]
