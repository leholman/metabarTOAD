#' PrimerStrip
#'
#'This function uses cutadapt to strip off non-biological nucleotides (adapters or primers) from fastq files. It depends on Cutadapt.
#'
#'@param PrimerF The 5-3 sequence of the forward primer. Ambiguous sequences are accepted (eg.N).
#'@param PrimerR The 5-3 sequence of the reverse primer. Ambiguous sequences are accepted (eg.N).
#'@param MinLen Minimum sequence length to be retained.
#'@param MaxLen Maximum sequence length to be retained.
#'@param UsePrimerFile TRUE or FALSE, is a csv file provided containing primer information?
#'@param folderwfiles A directory indicating files to be processed if not default folder.
#'@param folderoutput Output directory.
#'@param cutadaptdest Full directory for the cutadapt executable if not in PATH.
#'@param ncores Number of processsor cores to use.
#'
#'
#'@return None
#'
#'@importFrom utils read.csv
#'@importFrom seqinr c2s comp s2c
#'
#'@export


PrimerStrip <- function(PrimerF=NA,
                        PrimerR=NA,
                        MinLen=NA,
                        MaxLen=NA,
                        UsePrimerFile=FALSE,
                        folderwfiles="2.mergedreads",
                        folderoutput="3.strippedreads",
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
    message("Using primer data and metadata file to trim primers and length truncate.")
  }
  if(!UsePrimerFile){
    ##Write some further checks in here for the primers
    if (!any(c(is.numeric(MinLen),is.numeric(MaxLen)))){stop("Min and max length should be numeric, check the inputs.")}
    sampleindex <- gsub("(^.*).merged.fastq","\\1",list.files(folderwfiles,pattern="*.fastq"))
    primerindex <- rep("dummyprimer",length(sampleindex))
    primerdata <- data.frame("PrimerPair" = "dummyprimer","F" = PrimerF,"R" = PrimerR,"minsize" = MinLen,"maxsize" = MaxLen)
    primers <- "dummyprimer"
  }

    for (primer in primers){
      loopforward <- gsub("I","N",primerdata$F[primerdata$PrimerPair==primer])
      loopreverse <- gsub("I","N",primerdata$R[primerdata$PrimerPair==primer])
      loopreverse.revcomp <- c2s(rev(comp(s2c(loopreverse),ambiguous = TRUE)))
      loopmin <- primerdata$minsize[primerdata$PrimerPair==primer]
      loopmax <- primerdata$maxsize[primerdata$PrimerPair==primer]

      for (loopsample in sampleindex[primerindex==primer]){
        #this argument needs work
        message(paste0("Merging ",loopsample))
        cutadaptarg <- paste("-g  ^",loopforward,"...",loopreverse.revcomp," -m ",loopmin," -M ",loopmax," -n 2 -j ",ncores," --discard-untrimmed -o ",folderoutput,"/",loopsample,".stripped.fastq ",folderwfiles,"/",loopsample,".merged.fastq",sep="")
        log <- system2(cutadaptdest,cutadaptarg,stdout = TRUE,stderr = TRUE)
        cat(file="log.txt", log , append=T, sep="\n")
      }
    }
  }




