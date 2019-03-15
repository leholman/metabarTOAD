#' RenameFiles
#'
#'This function renames .fq files according to a user supplied index.
#'
#'@return None
#'
#'
#'@export


RenameFiles <- function(SeqIDs,DesiredIDs,folderwfiles="1.rawreads",usemetadatafile=FALSE){
  if(usemetadatafile==TRUE){
    index <- read.csv("metadata.csv")} else {
  if(length(SeqIDs)!=length(DesiredIDs)){stop("Sequencing and Desired IDs vary in length, check the inputs.")
      } else {index <- data.frame(RunID=SeqIDs,RealID=DesiredIDs)
  }}
  files <- list.files("1.rawreads/",pattern=".fastq")
  for (i in 1:length(files)){
    currentfile <- files[i]
    if(length(strsplit(currentfile,"_")[[1]])>5){stop("File name in unsupported format, check the sample idntifier does not contain underscores.")}
    currentfileID <- strsplit(currentfile,"_")[[1]][1]
    if(!(currentfileID %in% as.character(index$RunID))){stop(paste("Sample",currentfileID,"not found in supplied IDs"))}
    paste(" 1.rawreads/",as.character(index$RealID[index$RunID==currentfileID]),"_",paste(strsplit(currentfile,"_")[[1]][2:5],collapse ="_"),sep="")
    command <-paste("1.rawreads/",currentfile,paste(" 1.rawreads/",as.character(index$RealID[index$RunID==currentfileID]),"_",paste(strsplit(currentfile,"_")[[1]][2:5],collapse ="_"),sep=""),sep="")
    system2('mv',args=command, stdout = NULL,stderr = NULL)
  }
}
