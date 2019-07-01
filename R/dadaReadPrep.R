#'dadaReadPrep
#'
#'This function prepares paired-end Illumina data for the DADA2 pipeline.
#'It takes in paired-end demultiplexed amplicon fastq data and primer data and uses Cutadapt to trim off the primer region.
#'Input files are expected to have the structure SampleID_S*_L001_R1_001.fastq.
#'
#'@param PrimerF Sequence of the forward primer
#'@param PrimerR Sequence of the reverse primer
#'@param UsePrimerFile TRUE or FALSE, is a csv file provided containing primer information?
#'@param folderwfiles A directory indicating files to be processed if not default folder.
#'@param folderoutput Output directory.
#'@param lookforsecondprimer TRUE or FALSE, should Cutadapt look for a linked primer in the sequence, this is useful if the single read length covers the entire amplicon.
#'@param cutadaptdest Full directory for the cutadapt executable if not in PATH.
#'@param ncores Number of processsor cores to use.
#'@return None
#'
#'@importFrom utils read.csv
#'@importFrom seqinr c2s comp s2c
#'
#'@export


dadaReadPrep <- function(PrimerF=NA,
                         PrimerR=NA,
                         UsePrimerFile=FALSE,
                         folderwfiles="1.rawreads",
                         folderoutput="7.DADA2/trimmed.reads",
                         lookforsecondprimer=FALSE,
                         cutadaptdest="cutadapt",
                         ncores=1){

  if(UsePrimerFile){
    if(!file.exists("metadata.csv")){stop("No metadata file found, metadata needed to link primers and samples.")}
    if(!file.exists("primers.csv")){stop("No primer data found.")}
    primerdata <- read.csv("primers.csv")
    primers <- as.character(primerdata$PrimerPair)
    index <- read.csv("metadata.csv")
    primerindex <- as.character(index$Primer)
    sampleindex <- as.character(index$RealID)
    if(length(list.files(folderwfiles,pattern="L001_R1_001.fastq"))!=length(list.files(folderwfiles,pattern="L001_R2_001.fastq"))){stop("Different number of forward and reverse reads, function aborted.")}
    #check the below two lines they feel dicey
    fileindex <- gsub("(^.*)_L001_R1_001.fastq(.gz)?","\\1",list.files(folderwfiles,pattern="L001_R1_001.fastq"))
    fileindex <-  fileindex[match(sampleindex,sapply(strsplit(basename(fileindex), "_"), `[`, 1))]
    if (!dir.exists("7.DADA2/trimmed.reads")){dir.create("7.DADA2/trimmed.reads")}
    message("Using primer data and metadata file to trim primers")
  }
  if(!UsePrimerFile){
    message("Using function supplied primer data to trim primers")
    ##Write some further checks in here for the primers
    fileindex <- gsub("(^.*)_L001_R1_001.fastq(.gz)?","\\1",list.files(folderwfiles,pattern="L001_R1_001.fastq"))
    sampleindex <- sapply(strsplit(basename(fileindex), "_"), `[`, 1)
    primerindex <- rep("dummyprimer",length(sampleindex))
    primerdata <- data.frame("PrimerPair" = "dummyprimer","F" = PrimerF,"R" = PrimerR)
    primers <- "dummyprimer"
    if (!dir.exists("7.DADA2/trimmed.reads")){dir.create("7.DADA2/trimmed.reads")}
  }

  for (primer in primers){
    loopforward <- gsub("I","N",primerdata$F[primerdata$PrimerPair==primer])
    loopreverse <- gsub("I","N",primerdata$R[primerdata$PrimerPair==primer])
    loopforward.revcomp <- c2s(rev(comp(s2c(loopforward),ambiguous = TRUE)))
    loopreverse.revcomp <- c2s(rev(comp(s2c(loopreverse),ambiguous = TRUE)))
    count <- 1

    if(lookforsecondprimer){
      for (loopsample in sampleindex[primerindex==primer]){
        message(paste0("Trimming primers from sample ",count," of ",length(sampleindex[primerindex==primer])," Sample Name: ",loopsample))
        forreadarg <- paste0("^",loopforward,"...",loopreverse.revcomp)
        revreadarg <- paste0("^",loopreverse,"...",loopforward.revcomp)
        cutadaptarg <- paste0("-a  ",forreadarg," -A ",revreadarg," -j ",ncores,
                           " --discard-untrimmed -o ",folderoutput,"/",loopsample,"_R1_stripped.fastq.gz -p ",
                           folderoutput,"/",loopsample,"_R2_stripped.fastq.gz ",getwd(),"/",folderwfiles,"/",
                           fileindex[sampleindex==loopsample],"_L001_R1_001.fastq.gz ",getwd(),"/",folderwfiles,"/",
                           fileindex[sampleindex==loopsample],"_L001_R2_001.fastq.gz")
        log <- system2(cutadaptdest,cutadaptarg,stdout = TRUE,stderr = TRUE)
        cat(file="log.txt", log , append=T, sep="\n")
        count <- count + 1
      }
    }
    if(!lookforsecondprimer){
      for (loopsample in sampleindex[primerindex==primer]){
        message(paste0("Trimming primers from sample ",count," of ",length(sampleindex[primerindex==primer])," Sample Name:",loopsample))
        cutadaptarg <- paste0("-g  ^",loopforward," -G ^",loopreverse," -j ",ncores,
                              " --discard-untrimmed -o ",folderoutput,"/",loopsample,"_R1_stripped.fastq.gz -p ",
                              folderoutput,"/",loopsample,"_R2_stripped.fastq.gz ",getwd(),"/",folderwfiles,"/",
                              fileindex[sampleindex==loopsample],"_L001_R1_001.fastq.gz ",getwd(),"/",folderwfiles,"/",
                              fileindex[sampleindex==loopsample],"_L001_R2_001.fastq.gz")
        log <- system2(cutadaptdest,cutadaptarg,stdout = TRUE,stderr = TRUE)
        cat(file="log.txt", log , append=T, sep="\n")
        count <- count + 1
      }
    }
  }

}
