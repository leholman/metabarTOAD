#'PoolNFilterReads
#'
#'This function pools and performs quality filtering on merged metabarcoding reads.
#'
#'@param FastqMaxEE Maximum expected error of pooled sequences
#'@param folderwfiles Folder containing fastq files to be processed. Default is 3.strippedreads.
#'@param folderoutput Output folder for processed reads. Default is 4.pooledsamples.
#'@param UsePrimerFile TRUE or FALSE, is a csv file provided containing primer information?
#'@param vsearchdest Specify the location of the vsearch executable if PATH not configured appropriately.
#'
#'
#'@return None
#'
#'
#'@export


PoolNFilterReads <- function(FastqMaxEE=1,
                             folderwfiles="3.strippedreads",
                             folderoutput="4.pooledsamples",
                             UsePrimerFile=FALSE,
                             vsearchdest="vsearch"){
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
    message(paste("Pooling",length(list.files(folderwfiles,pattern=".fastq")),"fastq files"))
    sampleindex <- gsub("(^.*).stripped.fastq","\\1",list.files(folderwfiles,pattern="*.fastq"))
    primerindex <- rep("AllSamples",length(sampleindex))
    primers <- "AllSamples"
  }

  for (primer in primers){
   if(!primer=="AllSamples"){message(paste0("Processing ",primer," primers"))}
    message("Pooling samples")
    files <- paste(folderwfiles,"/",sampleindex[primerindex==primer],".stripped.fastq",sep="")
    catoutput <- paste(paste(files, collapse = " ")," > ",folderoutput,"/",primer,".pooled.fastq",sep="")
    system2("cat",catoutput)
    message("Done")
    qualityarg <- paste("-fastq_filter",paste(folderoutput,"/",primer,".pooled.fastq",sep=""),"--fastq_qmax 42 -fastq_maxee",FastqMaxEE," -fastaout",paste(folderoutput,"/",primer,".pooled.QF.fastq",sep=""),sep=" ")
    message("Filtering out poor quality sequences")
    log <- system2(vsearchdest, args=qualityarg,stdout = TRUE,stderr = TRUE)
    cat(file="log.txt", log , append=T, sep="\n")
    singletonsarg <- paste("-derep_fulllength",paste(folderoutput,"/",primer,".pooled.QF.fastq",sep=""),"-sizeout -output",paste(folderoutput,"/",primer,".pooled.ST.QF.fastq",sep=""),sep=" ")
    message("Done")
    message("Filtering out singletons")
    log <- system2(vsearchdest, args=singletonsarg,stdout = TRUE,stderr = TRUE)
    cat(file="log.txt", log , append=T, sep="\n")
    message("Done")
  }
}


