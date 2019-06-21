#'DenoiseUNOISE3
#'
#'This function denoises amplicons using the UNOISE3 algorithm from the USEARCH softwae by Robert Edgar.
#'It depends on the folder structure created by the function Folders().
#'
#'
#'@param usearchdest Specify the location of the usearch executable if PATH not configured appropriately.
#'@param UsePrimerFile TRUE or FALSE, is a csv file provided containing primer information?
#'
#'@return None
#'
#'@importFrom utils read.csv write.table
#'@importFrom seqinr read.fasta
#'
#'@export


DenoiseUNOISE3 <- function(UsePrimerFile=FALSE,usearchdest="usearch"){

  if(UsePrimerFile){
    if(!file.exists("metadata.csv")){stop("No metadata file found, metadata needed to link primers and samples.")}
    if(!file.exists("primers.csv")){stop("No primer data found.")}
    primerdata <- read.csv("primers.csv")
    primers <- as.character(primerdata$PrimerPair)
    index <- read.csv("metadata.csv")
    primerindex <- as.character(index$Primer)
    sampleindex <- as.character(index$RealID)
    message("Using primer data and metadata file.")
  }
  if(!UsePrimerFile){
    sampleindex <- gsub("(^.*).stripped.fastq","\\1",list.files("3.strippedreads",pattern="*.fastq"))
    primerindex <- rep("AllSamples",length(sampleindex))
    primers <- "AllSamples"
  }
  for (primer in primers){
    message(paste("Denoising OTUs"))
    unoisearg <- paste("-unoise3 4.pooledsamples/",primer,".pooled.ST.QF.fastq -zotus 5.OTUs/",primer,".unoise3.OTUs.fasta",sep="")
    log <- system2(usearchdest,unoisearg,stdout = TRUE,stderr = TRUE)
    cat(file="log.txt", log , append=T, sep="\n")
    loopfiles <- sampleindex[primerindex==primer]
    OTUs <- read.fasta(paste("5.OTUs/",primer,".unoise3.OTUs.fasta",sep=""))
    results <-data.frame("OTU"=paste("Zotu",1:length(OTUs),sep=""))
    blank <- data.frame(matrix(rep("0",length(OTUs)*length(loopfiles)),nrow = length(OTUs), ncol=length(loopfiles)))
    results <-cbind(results,blank)
    #this expression turns all the columns except the first into numeric variables
    results <- cbind(results[,1],data.frame(lapply(results[,2:(length(loopfiles)+1)], function(x) as.numeric(as.character(x)))))
    names(results) <- c("OTU",loopfiles)
    for (loopfile in loopfiles){
      message(paste("Mapping",loopfile,"to denoised seqs"))
      if(file.info(paste0("3.strippedreads/",loopfile,".stripped.fastq"))$size == 0){next}
      maparg <- paste("-usearch_global 3.strippedreads/",loopfile,".stripped.fastq -db 5.OTUs/",primer,".unoise3.OTUs.fasta -id 0.97 -maxaccepts 8 -maxrejects 256 -blast6out 6.mappings/hits.unoise3.",primer,".",loopfile,".txt -strand plus -maxhits 1",sep="")
      log <- system2(usearchdest,maparg,stdout = TRUE,stderr = TRUE)
      cat(file="log.txt", log , append=T, sep="\n")
      if(file.size(paste("6.mappings/hits.unoise3.",primer,".",loopfile,".txt",sep="")) == 0){next}else{
        hits <- read.csv(paste("6.mappings/hits.unoise3.",primer,".",loopfile,".txt",sep=""),sep="\t", header=F, stringsAsFactors=F)
        names(hits) <- c("query", "otu", "ident", "length", "mism", "gap", "qstart", "qend", "target_s", "target_e", "e.value", "bitscore")
        hitstab <- data.frame(table(hits$otu))
        for (item in 1:length(hitstab$Var1)){
          results[match(hitstab$Var1[item],results$OTU),loopfile] <- hitstab$Freq[item]
        }
      }
    }
    #Write out raw OTU x sample table UNOISE3
    message("mapping finished")
    rownames(results) <- results[,1]
    results <- results[,-1]
    write.table(results,file=paste("6.mappings/OTUtabs/",primer,".raw.unoise3.csv",sep=""),sep=",")
    message("OTU by sample table written.")
    message("Done")
  }
}






