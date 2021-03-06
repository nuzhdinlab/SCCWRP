library("plyr")
library(dplyr)
library("ggplot2")
library(lubridate)
library("ape")
library("vegan")
library("microbiome")
library(data.table)
library(tidyr)
library(hierDiversity)

#This script focuses on generating co-occurrence networks on a HUC-8 watershed scale within the SCCWRP archive.
#Subsetting waterhedss by land use bands to check for uniformity of co-occurrence network formation
#within similar watersheds to changes in land use.
setwd("~/Desktop/SCCWRP")
#Read in site data containing biological counts, water chemistry, and land usage
#values.  If this file is not yet generated then proceed with the following commands
#to generate it in the first place.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Filter out to taxonomic groups of interest.
GISBioData <- subset(GISBioData, MeasurementType == "benthic macroinvertebrate relative abundance")
#Read in sample metadata.
SCCWRP <- read.table("CSCI.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Get samples per watershed.
watersheds <- as.data.frame((table(SCCWRP$Watershed)))
colnames(watersheds) <- c("Watershed","Samples")
#Get the samples per watershed for watersheds with at least a certain number of samples.
LargeWatersheds <- subset(watersheds,Samples>=40)

#Initialize an empty gamma diversity data frame.
gammaDiversity <- data.frame()

for(WS in LargeWatersheds$Watershed){
  #Get all samples for a given watershed.
  GISBioDataWS <- subset(GISBioData,Watershed==WS)
  #print(paste(WS,median(GISBioDataWS$LU_2000_5K),median(GISBioDataWS$Ag_2000_5K),median(GISBioDataWS$CODE_21_2000_5K),median(GISBioDataWS$URBAN_2000_5K)))
  #Set land use bounds.
  LUmin <- 0
  LUmax <- 100
  #Subset samples per watershed within a land use band.
  GISBioDataWSLU <- subset(GISBioDataWS,LU_2000_5K >= LUmin & LU_2000_5K < LUmax)
  #Order samples by land use.
  GISBioDataWSLU <- arrange(GISBioDataWSLU,LU_2000_5K)
  #Get the number of samples which fit the land use profile for a given watersheds.
  nSamples <- length(unique(GISBioDataWSLU$UniqueID))
  nMin <- 20
  #Get the samples with the lowest land use value within the watershed by land use bands.
  GISBioDateWSLU <- GISBioDataWSLU[GISBioDataWSLU$UniqueID %in% as.vector(unique(GISBioDataWSLU$UniqueID)[1:nMin]),]
  if(nSamples >= nMin) {
    meanLU <- mean(na.omit(GISBioDataWSLU$LU_2000_5K))
    #Initialize a data frame where the rows are all of the unique measurements for a given
    #subset of the data.
    #Order the data frame by measurement name.
    selected <- arrange(GISBioDataWSLU,Year,UniqueID)
    eLSAInput <- as.data.frame(unique(selected$FinalID))
    colnames(eLSAInput)<-c("FinalID")
    eLSAInput <- as.data.frame(eLSAInput[order(as.character(eLSAInput$FinalID)),])
    colnames(eLSAInput)<-c("FinalID")
    #Plot CDF of land use for a watershed.
    #plot(ecdf(selected$LU_2000_5K),xlim=c(0.0,100.0),xlab="Land Use %cover",ylab="CDF",main="CDF for Land Use")
    #par(new=T)
    print(paste(WS,summary(ecdf(selected$LU_2000_5K))[3:4]))
    
    #Add the relative taxa abundances by column to a new dataframe.
    #The rows are the unique taxa in a given subset of data.
    selected <- selected[order(selected$Year,selected$UniqueID,selected$FinalID),]
    for(ID in unique(selected$UniqueID)){
      tmp <- filter(selected, UniqueID == ID)[,c("FinalID","Measurement","UniqueID")]
      tmp <- as.data.frame(tmp[order(tmp$FinalID),])
      tmp <- tmp[-c(3)]
      colnames(tmp)<-c("FinalID",paste("Measurement",ID,sep=" "))
      eLSAInput <- join(eLSAInput,tmp,by="FinalID")
      eLSAInput$FinalID=as.character(eLSAInput$FinalID)
      eLSAInput <- eLSAInput %>% group_by(FinalID) %>% summarise_if(is.numeric,mean,na.rm=TRUE)
      #print(ID)
    }
    
    eLSAInput[is.na(eLSAInput)] <- "NA"
    
    #Determine the number of time points in the eLSA input file.
    spotNum = length(unique(selected$Year))
    #Determine the number of replicates per time point in the eLSA input file.
    #In order to ensure a uniform number of replicates per year this needs to
    #be the maximum number of replicates for all of the years available.
    repMax = 0
    for(year in unique(selected$Year)){
      tmp <- filter(selected, Year == year)[,c("UniqueID","Year")]
      repNum = length(unique(tmp$UniqueID))
      if(repNum >= repMax){repMax = repNum}
      #print (paste(repMax,repNum,year,sep=" "))
    }
    repNum = repMax
    
    #Now insert the replicates with actual data in between the "NA" dummy columns
    #which ensure that the final eLSA input file has an even number of replicates
    #per year regardless of the variations in the actual number of sites (replicates)
    #sampled per year.
    eLSAtmp <- eLSAInput[,1]
    j=1
    k=1
    nulCol <- data.frame(matrix(ncol=repNum*spotNum,nrow=length(unique(selected$FinalID))))
    nulCol[,1] <- eLSAInput[,1]
    for(year in unique(selected$Year)){
      tmp <- filter(selected, Year == year)
      rep = length(unique(tmp$UniqueID))
      for(i in 1:repNum){
        if(i <= rep){
          repLabel = paste(year,"DoneRep",i,sep="")
          j=j+1
          k=k+1
          eLSAtmp[,k] <- eLSAInput[,j]
        }
        else{
          repLabel = as.character(paste(year,"Rep",i,sep=""))
          k=k+1
          eLSAtmp[,k] <- "NA"
          #print(paste(k,repLabel,sep=" "))
        }
      }
    }
    
    eLSAInput <- eLSAtmp
    #Designate a unique filename.
    #N is the number of samples in the subsample group.
    #S is the number of spots, or years represented in the subsample group.
    #R is the number of replicates per year.  Many of the years will have null replicates, but a uniform number is needed for eLSA.
    #M is the mean LU_2000_5K score per subsample group.
    filename = paste("SCCWRPWSbyLUWatershed",gsub(" ","",WS,fixed=TRUE),nMin,"Samples","LU",LUmin,"to",LUmax,"S",spotNum,"R",repNum,"M",meanLU,sep="")
    #Output file for use in eLSA.
    #write.table(eLSAInput,paste(filename,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    eLSACommand = paste("lsa_compute ",filename,".txt ","-r ",repNum," -s ",spotNum," ",filename,"Network.txt;",sep="")
    #print(eLSACommand)
    #Create a community matrix to determine gamma diversity.
    #Gamma0 = sample group richness, Gamma1 = sample group Shannon index, Gamma2 = sample group inverse Simpson index.
    abundances <- eLSAtmp[,-c(1)]
    abundances[] <- lapply(abundances,gsub,pattern="NA",replacement=as.numeric(0),fixed=TRUE)
    abundances <- as.data.frame(sapply(abundances,as.numeric))
    abundances <- t(abundances)
    gamma0 <- dz(abundances,lev="gamma",q=0)
    gamma1 <- dz(abundances,lev="gamma",q=1)
    gamma2 <- dz(abundances,lev="gamma",q=2)
    #print(paste(filename,gamma0,gamma1,gamma2))
    dat <- data.frame()
    dat[1,1] <- paste(filename,"Network.txt",sep="")
    dat[1,2] <- gamma0
    dat[1,3] <- gamma1
    dat[1,4] <- gamma2
    gammaDiversity <- rbind(gammaDiversity,dat)
  }
}

colnames(gammaDiversity) <- c("filename","gamma0","gamma1","gamma2")

#Read in eLSA output.
#Compute network statistics of the likeliest association networks between taxa.
library(igraph)
library(network)
library(stringr)
#Read in site data.
GISBioData <- read.table("CAGISBioData.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Ensure that all sites have a land use value.
GISBioData <- subset(GISBioData, LU_2000_5K != "NA")
#Filter out to taxonomic groups of interest.
GISBioData <- subset(GISBioData, MeasurementType == "benthic macroinvertebrate relative abundance")
networkfiles <- Sys.glob("SCCWRPWSbyLUWatershed*20SamplesLU*Network.txt")
networkAnalysis <- data.frame()
networkConTaxa <- data.frame()
networkCovTaxa <- data.frame()
#Define a 'not in' function.
'%!in%' <- function(x,y)!('%in%'(x,y))
for(networkFile in networkfiles){
  networkdata <- read.delim(networkFile,header=TRUE, sep="\t",as.is=T,check.names=FALSE)
  #Filter out association network data based on P and Q scores, for the local similarity
  #between two factors, with values less than a particuar threshold.
  networkdata <- filter(networkdata, P <= 1e-2)
  networkdata <- filter(networkdata, Q <= 1e-2)
  names(networkdata)[names(networkdata)=="LS"]<-"weight"
  meanLU <- as.numeric(str_match(networkFile, "SCCWRPWSbyLUWatershed(.*?)20SamplesLU(.*?)to(.*?)S(.*?)R(.*?)M(.*?)Network.txt")[7])
  #Generate network graph and begin calculating network parameters.
  networkgraph=graph.data.frame(networkdata,directed=FALSE)
  Network_size<-network.size(as.network(get.adjacency(networkgraph,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency"))
  if(ecount(networkgraph)>0){
    #Get the full weighted adjacency matrix.
    networkmatrix <- as.matrix(get.adjacency(networkgraph,attr='weight'))
    #Mean interaction strength
    meanStrength <- mean(abs(networkmatrix))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_network <- eigen(networkmatrix)
    #Get the real component first eigenvalue.
    lambda_network_m <- Re(lambda_network$values[1])
    #Generate randomized version of full weighted adjacency matrix.
    set.seed(1)
    randnetworkmatrix <- matrix(sample(as.vector((networkmatrix))),nrow=nrow(networkmatrix),ncol=ncol(networkmatrix))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_rand <- eigen(randnetworkmatrix)
    #Get the real component of the first eigenvalue.
    lambda_rand_m <- Re(lambda_rand$values[1])
    #Calculate stability parameter.
    gamma <- lambda_network_m/lambda_rand_m
    #Calculate the degree heterogeneity.
    networkmatrix[upper.tri(networkmatrix)] <- 0
    networkmatrix <- ifelse(networkmatrix!=0,1,networkmatrix)
    zeta <- mean(colSums(networkmatrix)^2)/mean(colSums(networkmatrix))^2
    #Calculate modularity
    networkModularity <- modularity(cluster_edge_betweenness(networkgraph, weights=NULL,directed=FALSE))
    M <- networkModularity
    networkNodecount <-network.size(as.network(get.adjacency(networkgraph,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency"))
    # Get the number of unique network edges
    networkEdgecount <- network.edgecount(as.network(get.adjacency(networkgraph,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency"))
    # Get the number of nodes
    networkNodecount <- network.size(as.network(get.adjacency(networkgraph,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency"))
    # Get the average degree per node.
    k <- (2*networkEdgecount)/networkNodecount
    # Calculate the modularity of the random network.
    networkRandModularity <- (1-(2/sqrt(networkNodecount)))*(2/k)^(2/3)
    # Calculate the log ratio of the modularities.
    l_rM <- log(networkModularity/networkRandModularity)
  }
  #Filter contravariant network data based on local similarity scores.
  networkdataCon <- subset(networkdata,networkdata$weight<0)
  #Aggregate significantly contravarying taxa.
  networkdataConTemp <- networkdataCon[,c("X","Y","weight")]
  networkdataConTemp <- as.data.frame(table(append(networkdataConTemp$X,networkdataConTemp$Y,after=length(networkdataConTemp$X))))
  networkdataConTemp$meanLU <- meanLU
  networkConTaxa <- rbind(networkConTaxa,networkdataConTemp)
  #Generate network graph and begin calculating network parameters.
  networkgraphCon=graph.data.frame(networkdataCon,directed=FALSE)
  if(ecount(networkgraphCon)>0){
    #Get the full weighted adjacency matrix.
    networkmatrix <- as.matrix(get.adjacency(networkgraphCon,attr='weight'))
    #Mean interaction strength
    meanStrength_Con <- mean(abs(networkmatrix))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_network <- eigen(networkmatrix)
    #Get the real component first eigenvalue.
    lambda_network_m_Con <- Re(lambda_network$values[1])
    #Generate randomized version of full weighted adjacency matrix.
    set.seed(1)
    randnetworkmatrix <- matrix(sample(as.vector((networkmatrix))),nrow=nrow(networkmatrix),ncol=ncol(networkmatrix))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_rand_Con <- eigen(randnetworkmatrix)
    #Get the real component of the first eigenvalue.
    lambda_rand_Con <- Re(lambda_rand_Con$values[1])
    #Calculate stability parameter.
    gamma_Con <- lambda_network_m_Con/lambda_rand_Con
    #Calculate the degree heterogeneity.
    networkmatrixCon <- networkmatrix
    networkmatrixCon[upper.tri(networkmatrixCon)] <- 0
    networkmatrixCon <- ifelse(networkmatrixCon!=0,1,networkmatrixCon)
    zeta_Con <- mean(colSums(networkmatrixCon)^2)/mean(colSums(networkmatrixCon))^2
    #Calculate the degree heterogeneity of the corresponding random network.
    randnetworkmatrixCon <- randnetworkmatrix
    randnetworkmatrixCon[upper.tri(randnetworkmatrixCon)] <- 0
    randnetworkmatrixCon <- ifelse(randnetworkmatrixCon!=0,1,randnetworkmatrixCon)
    zeta_rand_Con <- mean(colSums(randnetworkmatrixCon)^2)/mean(colSums(randnetworkmatrixCon))^2
    # Log response ratio of degree heterogeneity.
    l_con_rzeta <- log(zeta_Con/zeta_rand_Con)
    # Generate adjacency matrix of relative taxa abundance correlations
    adj= as.network(get.adjacency(networkgraphCon,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency")
    # Get the number of unique network edges
    networkEdgecount <- network.edgecount(adj)
    networkEdgecountCon <- networkEdgecount
    # Get the number of nodes
    networkNodecount <- network.size(adj)
    # Get the average degree per node.
    k <- (2*networkEdgecount)/networkNodecount
    # Get the random characteristic path length.
    networkRandLength <- 0.5+((log(networkNodecount)-0.5772156649)/log(k))
    # Get the random clustering coefficient.
    networkRandClustering <- k/networkNodecount
    # Get the network density.
    networkDensity <- network.density(adj)
    con_C <- networkDensity
    # Calculate the modularity of the network.
    networkModularity <- modularity(cluster_edge_betweenness(networkgraphCon, weights=NULL,directed=FALSE))
    con_M <- networkModularity
    # Calculate the number of groups related to the modularity value.
    networkModGroups <- length(cluster_edge_betweenness(networkgraphCon, weights=NULL,directed=FALSE))
    # Calculate the average network path length
    networkLength <- mean_distance(networkgraphCon,directed=FALSE)
    con_L <- networkLength
    # Calculate the clustering coefficient
    networkClustering <- transitivity(networkgraphCon,type="globalundirected",isolate="zero")
    con_Cl <- networkClustering
    # Calcuate the log ratio of clustering coefficients.
    l_con_rCl <- log(networkClustering/networkRandClustering)
    # Calculate the modularity of the random network.
    networkRandModularity <- (1-(2/sqrt(networkNodecount)))*(2/k)^(2/3)
    # Calculate the log ratio of the modularities.
    l_con_rM <- log(networkModularity/networkRandModularity)
    # Get log ratio of characteristic path lengths.
    l_con_rL <- log(networkLength/networkRandLength)
  }
  #Filter covariant network data based on local similarity scores.
  networkdataCov <- subset(networkdata,networkdata$weight>0)
  #Aggregate significantly contravarying taxa.
  networkdataCovTemp <- networkdataCov[,c("X","Y","weight")]
  networkdataCovTemp <- as.data.frame(table(append(networkdataCovTemp$X,networkdataCovTemp$Y,after=length(networkdataCovTemp$X))))
  networkdataCovTemp$meanLU <- meanLU
  networkCovTaxa <- rbind(networkCovTaxa,networkdataCovTemp)
  #Generate network graph and begin calculating network parameters.
  networkgraphCov=graph.data.frame(networkdataCov,directed=FALSE)
  if(ecount(networkgraph)>0){
    #Get the full weighted adjacency matrix.
    networkmatrix <- as.matrix(get.adjacency(networkgraphCov,attr='weight'))
    #Mean interaction strength
    meanStrength_Cov <- mean(abs(networkmatrix))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_network <- eigen(networkmatrix)
    #Get the real component first eigenvalue.
    lambda_network_m_Cov <- Re(lambda_network$values[1])
    #Generate randomized version of full weighted adjacency matrix.
    set.seed(1)
    randnetworkmatrix <- matrix(sample(as.vector((networkmatrix))),nrow=nrow(networkmatrix),ncol=ncol(networkmatrix))
    #Get the eigenvalues of the full weighted adjacency matrix.
    lambda_rand_Cov <- eigen(randnetworkmatrix)
    #Get the real component of the first eigenvalue.
    lambda_rand_Cov <- Re(lambda_rand_Cov$values[1])
    #Calculate stability parameter.
    gamma_Cov <- lambda_network_m_Cov/lambda_rand_Cov
    #Calculate the degree heterogeneity.
    networkmatrixCov <- networkmatrix
    networkmatrixCov[upper.tri(networkmatrixCov)] <- 0
    networkmatrixCov <- ifelse(networkmatrixCov!=0,1,networkmatrixCov)
    zeta_Cov <- mean(colSums(networkmatrixCov)^2)/mean(colSums(networkmatrixCov))^2
    #Calculate the degree heterogeneity of the corresponding random network.
    randnetworkmatrixCov <- randnetworkmatrix
    randnetworkmatrixCov[upper.tri(randnetworkmatrixCov)] <- 0
    randnetworkmatrixCov <- ifelse(randnetworkmatrixCov!=0,1,randnetworkmatrixCov)
    zeta_rand_Cov <- mean(colSums(randnetworkmatrixCov)^2)/mean(colSums(randnetworkmatrixCov))^2
    # Log response ratio of degree heterogeneity.
    l_cov_rzeta <- log(zeta_Cov/zeta_rand_Cov)
    # Generate adjacency matrix of relative taxa abundance correlations
    adj= as.network(get.adjacency(networkgraphCov,attr='weight',sparse=FALSE),directed=FALSE,loops=FALSE,matrix.type="adjacency")
    # Get the number of unique network edges
    networkEdgecount <- network.edgecount(adj)
    networkEdgecountCov <- networkEdgecount
    # Get the number of nodes
    networkNodecount <- network.size(adj)
    # Get the average degree per node.
    k <- (2*networkEdgecount)/networkNodecount
    # Get the random characteristic path length.
    networkRandLength <- 0.5+((log(networkNodecount)-0.5772156649)/log(k))
    # Get the random clustering coefficient.
    networkRandClustering <- k/networkNodecount
    # Get the network density.
    networkDensity <- network.density(adj)
    cov_C <- networkDensity
    # Calculate the modularity of the network.
    networkModularity <- modularity(cluster_edge_betweenness(networkgraphCov, weights=NULL,directed=FALSE))
    cov_M <- networkModularity
    # Calculate the number of groups related to the modularity value.
    networkModGroups <- length(cluster_edge_betweenness(networkgraphCov, weights=NULL,directed=FALSE))
    # Calculate the average network path length
    networkLength <- mean_distance(networkgraphCov,directed=FALSE)
    cov_L <- networkLength
    # Calculate the clustering coefficient
    networkClustering <- transitivity(networkgraphCov,type="globalundirected",isolate="zero")
    cov_Cl <- networkClustering
    # Calcuate the log ratio of clustering coefficients.
    l_cov_rCl <- log(networkClustering/networkRandClustering)
    # Calculate the modularity of the random network.
    networkRandModularity <- (1-(2/sqrt(networkNodecount)))*(2/k)^(2/3)
    # Calculate the log ratio of the modularities.
    l_cov_rM <- log(networkModularity/networkRandModularity)
    # Get log ratio of characteristic path lengths.
    l_cov_rL <- log(networkLength/networkRandLength)
  }
  dat <- data.frame()
  dat[1,1] <- networkFile
  dat[1,2] <- meanLU
  dat[1,3] <- l_con_rL
  dat[1,4] <- l_con_rCl
  dat[1,5] <- l_con_rM
  dat[1,6] <- l_cov_rL
  dat[1,7] <- l_cov_rCl
  dat[1,8] <- l_cov_rM
  dat[1,9] <- lambda_network_m
  dat[1,10] <- con_L
  dat[1,11] <- con_Cl
  dat[1,12] <- con_M
  dat[1,13] <- cov_L
  dat[1,14] <- cov_Cl
  dat[1,15] <- cov_M
  dat[1,16] <- zeta
  dat[1,17] <- con_C
  dat[1,18] <- cov_C
  dat[1,19] <- Network_size
  dat[1,20] <- Pm <- networkEdgecountCov/(networkEdgecountCov+networkEdgecountCon)
  dat[1,21] <- lambda_network_m_Con
  dat[1,22] <- lambda_network_m_Cov
  dat[1,23] <- zeta_Con
  dat[1,24] <- zeta_Cov
  dat[1,25] <- gamma_Con
  dat[1,26] <- gamma_Cov
  dat[1,27] <- gamma
  dat[1,28] <- lambda_rand_m
  dat[1,29] <- lambda_rand_Con
  dat[1,30] <- lambda_rand_Cov
  dat[1,31] <- M
  dat[1,32] <- l_rM
  dat[1,33] <- meanStrength
  dat[1,34] <- meanStrength_Cov
  dat[1,35] <- meanStrength_Con
  dat[1,36] <- zeta_rand_Con
  dat[1,37] <- l_con_rzeta
  dat[1,38] <- zeta_rand_Cov
  dat[1,39] <- l_cov_rzeta
  dat[1,40] <- Watershed <- str_match(networkFile,"SCCWRPWSbyLUWatershed(.*?)20SamplesLU(.*?)to(.*?)S(.*?)R(.*?)M(.*?)Network.txt")[2]
  networkAnalysis <- rbind(networkAnalysis,dat)
  print(paste(networkFile,meanLU,l_con_rL,l_con_rCl,l_con_rM,l_cov_rL,l_cov_rCl,l_cov_rM,lambda_network_m,con_L,con_Cl,con_M,cov_L,cov_Cl,cov_M,zeta,con_C,cov_C,Network_size,Pm,lambda_network_m_Con,lambda_network_m_Cov,zeta_Con,zeta_Cov,gamma_Con,gamma_Cov,gamma,lambda_rand_m,lambda_rand_Con,lambda_rand_Cov,M,l_rM,meanStrength,meanStrength_Cov,meanStrength_Con,zeta_rand_Con,l_con_rzeta,zeta_rand_Cov,l_cov_rzeta,Watershed))
}
colnames(networkAnalysis) <- c("filename","meanLU","l_con_rL","l_con_rCl","l_con_rM","l_cov_rL","l_cov_rCl","l_cov_rM","lambda_network_m","con_L","con_Cl","con_M","cov_L","cov_Cl","cov_M","zeta","con_C","cov_C","Network_size","Pm","lambda_network_m_Con","lambda_network_m_Cov","zeta_Con","zeta_Cov","gamma_Con","gamma_Cov","gamma","lambda_rand_m","lambda_rand_Con","lambda_rand_Cov","M","l_rM","meanStrength","meanStrength_Cov","meanStrength_Con","zeta_rand_Con","l_con_rzeta","zeta_rand_Cov","l_cov_rzeta","Watershed")
networkAnalysis[networkAnalysis=="-Inf"] <- NA
networkAnalysis[networkAnalysis=="Inf"] <- NA
networkAnalysis <- arrange(networkAnalysis,meanLU)

write.table(networkAnalysis,"SCCWRPWatershedByLU20Samples.txt",quote=FALSE,sep="\t",row.names = FALSE)
networkAnalysis <- read.table("SCCWRPWatershedByLU20Samples.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

#Check how network topologies change within watersheds.
test <- networkAnalysis
test$diff_meanLU <- ave(test$meanLU, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_zeta_Cov <- ave(test$zeta_Cov, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_zeta_Con <- ave(test$zeta_Con, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_l_cov_rM <- ave(test$l_con_rM, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_l_con_rM <- ave(test$l_con_rM, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_cov_C <- ave(test$cov_C, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_con_C <- ave(test$con_C, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_meanStrength_Cov <- ave(test$meanStrength_Cov, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_meanStrength_Con <- ave(test$meanStrength_Con, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_lambda_network_m_Cov <- ave(test$lambda_network_m_Cov, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_lambda_network_m_Con <- ave(test$lambda_network_m_Con, test$Watershed, FUN=function(x) c(NA, diff(x)))
test$diff_lambda_network_m <- ave(test$lambda_network_m, test$Watershed, FUN=function(x) c(NA, diff(x)))
