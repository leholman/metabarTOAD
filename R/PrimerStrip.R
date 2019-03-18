#' PrimerStrip
#'
#'This function uses cutadapt to strip off non-biological nucleotides (adapters or primers) from fastq files. It depends on Cutadapt.
#'
#'@param PrimerF The 5-3 sequence of the forward primer. Ambiguous sequences are accepted (eg.N).
#'@param PrimerR The 5-3 sequence of the reverse primer. Ambiguous sequences are accepted (eg.N).
#'@return None
#'
#'@importFrom utils read.csv
#'
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
                        ncores=1){0
  if(UsePrimerFile){
    if(!file.exists("metadata.csv")){stop("No metadata file found, metadata needed to link primers and samples.")}
    if(!file.exists("primers.csv")){stop("No primer data found.")}
    PrimerData <- read.csv("primers.csv")
    Primers <- as.character(PrimerData$PrimerPair)
    index <- read.csv("metadata.csv")
    primerindex <- as.character(index$Primer)
    sampleindex <- as.character(index$RealID)
    message("Using primer data and metadata file to trim primers and length truncate.")
  }
  if(!UsePrimerFile){
    ##Write somr checks in here ot make sure the data handed to the fucntion is correct
    sampleindex <- gsub("(^.*).merged.fastq","\\1",list.files(folderwfiles,pattern="*.fastq"))
    primerindex <- rep("dummyprimer",length(sampleindex))
    PrimerData <- data.frame("PrimerPair" = "dummyprimer","F" = PrimerF,"R" = PrimerR,"minsize" = MinLen,"maxsize" = MaxLen)
    Primers <- "dummyprimer"
  }

    for (primer in Primers){
      loopforward <- gsub("I","N",PrimerData$F[PrimerData$PrimerPair==primer])
      loopreverse <- gsub("I","N",PrimerData$R[PrimerData$PrimerPair==primer])
      loopreverse.revcomp <- c2s(rev(comp(s2c(loopreverse),ambiguous = TRUE)))
      loopmin <- PrimerData$minsize[PrimerData$PrimerPair==primer]
      loopmax <- PrimerData$maxsize[PrimerData$PrimerPair==primer]

      for (sample in sampleindex[primerindex==primer]){
        #this argument needs work
        cutadaptarg <- paste("-g  ",loopforward,"...",loopreverse.revcomp," -m ",loopmin," -M ",loopmax," -n 2 -j ",ncores," --discard-untrimmed -o ",folderoutput,"/",sample,".stripped.fastq ",folderwfiles,"/",sample,".merged.fastq",sep="")

      }

    }


    #######YOU NEED TO FIX THIS PART
    ####HOW TO LOOP OVER MULTIPLE PRIMERS AND ALSO DO IT ONCE WHEN SET SIMPLY

  }



}
