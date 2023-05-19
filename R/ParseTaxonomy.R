#'ParseTaxonomy
#'
#'This function parses taxonomy from a BLAST query to provide assignments at a given identity.
#'An example blast command to use to generate a dataset to use for assignbments might be blastn -query ../path/to/your/seqs.fasta -db nt -out ../path/to/your/resultsfolder/results.txt -outfmt '6 qseqid qlen slen qcovs qcovhsp sseqid bitscore score evalue pident qstart qend sstart send staxids sscinames' -num_alignments 25 -num_threads 1
#'
#'
#'@param pctThreshold Percent identity match threashold for high quality assignments, 97-100 depending on marker
#'@param covpct Percent coverage of hit to be considered high quality
#'@param lwrpctThreshold Percent identity match threshold under which a hit is not considered
#'@param lwrcovpct Percent coverage of hit under which a hit is not considered
#'@param blastoutput The full directory to a text file with output from the BLAST standalone executable
#'@param lineages The full directory to a text file containing lineage information for each NCBI taxonomy number
#'@param returntable Return a table showing the assignments by catagory?
#'@return Returns a data frame object containing assignments
#'@importFrom data.table fread
#'@importFrom utils capture.output
#'
#'@export


ParseTaxonomy <- function(pctThreshold=97,
                          covpct=90,
                          lwrpctThreshold=85,
                          lwrcovpct=65,
                          blastoutput=NA,
                          lineages=NA,
                          returntable=TRUE){

  #Set up results table and read in data
  ## check files are real
  message("Reading in data")
  data <- fread(blastoutput,sep="\t",header=T)
  lineages <- fread(lineages,header=T)
  message("Done")
  #Data checks to go in here

  ## lineages most have OTU Kingdom -> Species
  ## All staxid must be in lineages
  ## pct Threashold should be 50-100 and bigger than lwrpctThreshold

  if (!all(names(data)==c("OTUID","qlen","slen","qcov","qcovhsp","sseqid","bitscore","score","evalue","pctid","qstart","qend","start","send","staxid","sscinames"))){
    data <- fread(blastoutput,sep="\t",header=F)
    names(data) <- c("OTUID","qlen","slen","qcov","qcovhsp","sseqid","bitscore","score","evalue","pctid","qstart","qend","start","send","staxid","sscinames")
  }
  assignmentResults <- data.frame("OTU"=as.character(unique(data$OTUID)),"assignmentQual"=rep(NA,length(as.character(unique(data$OTUID)))),"taxID"=rep(NA,length(as.character(unique(data$OTUID)))))


  #Step 1 clean raw data of unwanted assignments
  message("Parsing high quality assignments")
  data2 <- data[-grep("uncultured | Uncultured | environmental | environmental | construct | clone",data$sscinames),]

  ##Step 2 Subset for only high conf assignments
  data3 <- data2[-grep("sp\\.",data2$sscinames),]
  HCdata <- data3[data3$pctid > pctThreshold & data3$qcov > covpct,]


  #function to return the number of unique taxa hits per OTU within the high confidence limits
  numHits <- function(x){
    return(length(unique(HCdata$staxid[HCdata$OTUID==x])))
  }
  hits.dist <- lapply(unique(as.character(HCdata$OTUID)),numHits)
  hitmismatches <- data.frame("OTU"=unique(as.character(HCdata$OTUID)),"hits"=unlist(hits.dist))
  Hits1 <- HCdata[HCdata$OTUID %in% hitmismatches$OTU[hitmismatches$hits==1],]

  ##mark high confidence assignments
  for (OTU in unique(as.character(Hits1$OTUID))){
    assignmentResults$assignmentQual[assignmentResults$OTU==OTU] <- "High"
    assignmentResults$taxID[assignmentResults$OTU==OTU] <- as.character(unique(Hits1$staxid[Hits1$OTUID==OTU]))
  }

  ##Identify high confidence assingments with more than one match
  HQassigns <- as.character(assignmentResults$OTU[assignmentResults$assignmentQual=="High"])
  HQassigns <- HQassigns[!is.na(HQassigns)]

  HCdata2 <- data2[data2$pctid > pctThreshold & data2$qcov > covpct & data2$qcovhsp > covpct,]

  hits.dist2 <- lapply(unique(as.character(HCdata2$OTUID)),numHits)
  hitmismatches2 <- data.frame("OTU"=unique(as.character(HCdata2$OTUID)),"hits"=unlist(hits.dist2))
  Hits1.2 <- HCdata2[HCdata2$OTUID %in% hitmismatches2$OTU[hitmismatches2$hits>1],]
  Hits1.2 <- Hits1.2[!as.character(Hits1.2$OTUID) %in% HQassigns,]

  ##mark high confidence assignments with multiple hits
  for (OTU in unique(as.character(Hits1.2$OTUID))){
    assignmentResults$assignmentQual[assignmentResults$OTU==OTU] <- "High-MH"
  }

  #Reassign ambiguous taxonomic IDs (high)
  reassign <-as.character(assignmentResults$OTU[grep(";",assignmentResults$taxID)] )
  assignmentResults$assignmentQual[assignmentResults$OTU %in% reassign]<- "High-MH-S"
  message("Done")


  ##Subset data to retain poor but meaningful hits
  message("Parsing low confidence hits")
  poorData <- data2[data2$pctid > lwrpctThreshold & data2$qcov > lwrcovpct,]

  ##mark low confidence hits
  conditional <- unique(as.character(poorData$OTUID)) %in% as.character(assignmentResults$OTU[is.na(assignmentResults$assignmentQual)])
  for (OTU in unique(as.character(poorData$OTUID))[conditional]){
    assignmentResults$assignmentQual[assignmentResults$OTU==OTU] <- "Low"
    assignmentResults$taxID[assignmentResults$OTU==OTU] <- as.character(poorData$staxid[poorData$OTUID==OTU][match(max(poorData$bitscore[poorData$OTUID==OTU]),poorData$bitscore[poorData$OTUID==OTU])])
  }
  message("Done")
  #Reassign ambiguous taxonomic IDs (low)
  message("Marking mulitple hits for low confidence assignments")
  reassign <- as.character(assignmentResults$OTU[grep(";",assignmentResults$taxID)] )
  reassign <- reassign[!reassign %in% assignmentResults$OTU[assignmentResults$assignmentQual=="High-MH-S"]]
  assignmentResults$assignmentQual[assignmentResults$OTU %in% reassign]<- "Low-MH-S"
  message("Done")


  message("Setting up data for LCA assignments")
  # Get rid of NAs
  assignmentResults$assignmentQual[is.na(assignmentResults$assignmentQual)] <- "None"
  assignmentResults$taxID[is.na(assignmentResults$taxID)] <- "None"


  ## Use linages to pull out Kingdom -> Species for HQ and Low confidence hits
  #the below line was used in place of the if else ladder below in previous versions
  #taxIds <- c(unlist(strsplit(nums[grep(";",nums)],";")),nums[-grep(";",nums)],unique(Hits1.2$staxid[-grep(";",Hits1.2$staxid)]))

  nums <- assignmentResults$taxID
  nums <- nums[-grep("None",nums)]
  if (length(grep(";",nums))>0 & length(Hits1.2$staxid[grep(";",Hits1.2$staxid)])>0) {taxIds <- c(unlist(strsplit(nums[grep(";",nums)],";")),nums[-grep(";",nums)],unique(Hits1.2$staxid[-grep(";",Hits1.2$staxid)]))
  } else if (length(Hits1.2$staxid[grep(";",Hits1.2$staxid)])>0) {taxIds <- c(nums,unique(Hits1.2$staxid[-grep(";",Hits1.2$staxid)]))
  } else if (length(grep(";",nums))>0) {taxIds <- c(unlist(strsplit(nums[grep(";",nums)],";")),nums[-grep(";",nums)],unique(Hits1.2$staxid))
  } else {taxIds <- c(nums,unique(Hits1.2$staxid))}
  taxIds <- unique(taxIds)
  #taxIds <-c(as.numeric(assignmentResults$taxID[-grep(";",assignmentResults$taxID)]),unlist(strsplit(assignmentResults$taxID[grep(";",assignmentResults$taxID)],";")))
  #smolLineages <- lineages[lineages$tax_id %in% as.numeric(assignmentResults$taxID),]
  smolLineages <- lineages[lineages$tax_id %in% taxIds,]
  smolLineages2 <- smolLineages[,1:8]
  smolLineages3 <- as.data.frame(smolLineages2)
  colnames(smolLineages2)[-1]
  row.names(smolLineages3) <- smolLineages3$tax_id


  assignmentResults <- cbind(assignmentResults,matrix(ncol=7,nrow=length(assignmentResults[,1])))
  colnames(assignmentResults) <- c(colnames(assignmentResults)[1:3],colnames(smolLineages2)[2:8])

  message("Done")

  message("Assigning high confidence hits")
  #Assign HC hits
  for (OTU in as.character(assignmentResults$OTU[assignmentResults$assignmentQual=="High"])){
    runningID <- assignmentResults$taxID[assignmentResults$OTU==OTU]
    if(!runningID %in% smolLineages$tax_id){next}
    assignmentResults[assignmentResults$OTU==OTU,] <- c(as.matrix(assignmentResults[assignmentResults$OTU==OTU,1:3]),
                                                        as.matrix(smolLineages2[smolLineages2$tax_id==runningID,][,2:8]))

  }
  message("Done")
  message("Assigning low confidence hits")
  #Assign LC hits
  lowOTUs <- as.character(assignmentResults$OTU[match(as.character(assignmentResults$OTU[assignmentResults$assignmentQual=="Low"]),assignmentResults$OTU)])
  lowOTUs <- lowOTUs[assignmentResults$taxID[match(lowOTUs,as.character(assignmentResults$OTU))] %in% smolLineages$tax_id]

  for (OTU in lowOTUs){
    runningID <- assignmentResults$taxID[assignmentResults$OTU==OTU]
    assignmentResults[assignmentResults$OTU==OTU,] <- c(as.matrix(assignmentResults[assignmentResults$OTU==OTU,1:3]),
                                                        as.matrix(smolLineages2[smolLineages2$tax_id==runningID,][,2:8]))

  }
  message("Done")

  ## Determine lowest comon ancestor for multi hit assignments

  LCA <- function(taxa_to_check){
    ##INSERT checks on input 1) should be character or numerical 2) should be >1 number
    #run chekc for smolLineages
    runningdata <- cbind(data.frame("ID"=taxa_to_check,"Lineage_Avail"=rep("yes",length(taxa_to_check))),as.data.frame(matrix(nrow=length(taxa_to_check),ncol=7)))
    row.names(runningdata) <- runningdata[,1]
    runningdata$Lineage_Avail <- as.character(runningdata$Lineage_Avail)


    for (ID in as.character(runningdata[,1])){
      if (!ID %in% as.character(smolLineages$tax_id)){runningdata$Lineage_Avail[runningdata$ID==ID] <- "No"
      next}
      runningdata[ID,3:9] <- smolLineages3[ID,2:8]
    }
    runningdata <- runningdata[!is.na(runningdata$V1),]
    num.clades <- apply(runningdata[,3:9],2,function(x) length(unique(x)))

    #else if ladder returning LCA
    out <- c()
    if (num.clades[6]==1 & sum(num.clades[1:5])==5){out <- as.character(apply(runningdata[,3:9],2,function(x) unique(x))[1:6])
    } else if (num.clades[5]==1 & sum(num.clades[1:4])==4){out <-(as.character(apply(runningdata[,3:9],2,function(x) unique(x))[1:5]))
      } else if (num.clades[4]==1 & sum(num.clades[1:3])==3){ out <-(as.character(apply(runningdata[,3:9],2,function(x) unique(x))[1:4]))
        } else if (num.clades[3]==1 & sum(num.clades[1:2])==2){ out <- (as.character(apply(runningdata[,3:9],2,function(x) unique(x))[1:3]))
          } else if (num.clades[2]==1 & sum(num.clades[1])==1){ out <- (as.character(apply(runningdata[,3:9],2,function(x) unique(x))[1:2]))
            } else if (num.clades[1]==1){out <- as.character(apply(runningdata[,3:9],2,function(x) unique(x))[1])
    }
  ##Insert NAs to give an output of size 7
  if (length(out)==0){return(warning("Can't find taxa in database"))}
  output <- rep(NA,7)
  output[1:length(out)] <- out
  return(output)

  }

  message("Running LCA on Assignments")
  #HQ-MH
  for (OTU in as.character(assignmentResults$OTU[assignmentResults$assignmentQual=="High-MH"])){
    runningtaxa <- unique(HCdata2$staxid[HCdata2$OTUID==OTU])
    assignmentResults[assignmentResults$OTU==OTU,4:10] <- LCA(runningtaxa)
  }

  #HQ-MH-S
  for (OTU in as.character(assignmentResults$OTU[assignmentResults$assignmentQual=="High-MH-S"])){
    runningtaxa <- c(strsplit(assignmentResults$taxID[assignmentResults$OTU==OTU],";"))
    assignmentResults[assignmentResults$OTU==OTU,4:10] <- LCA(runningtaxa)
  }

  #Low-MH-S
  for (OTU in as.character(assignmentResults$OTU[assignmentResults$assignmentQual=="Low-MH-S"])){
    runningtaxa <- c(strsplit(assignmentResults$taxID[assignmentResults$OTU==OTU],";"))
    assignmentResults[assignmentResults$OTU==OTU,4:10] <- LCA(runningtaxa)
  }
  message("Done")

  ## clean up - check for no asignments and reallocate as none
  if(returntable==TRUE){message(paste0(capture.output(table(assignmentResults$assignmentQual)), collapse = "\n"))}

  return(assignmentResults)
}

