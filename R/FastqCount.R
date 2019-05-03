#' FastqCount
#'
#'This function counts the number of reads in a fastq file. It handles one file at a time.
#'
#'@param file The file location.
#'
#'@return The number of reads in the fastq file
#'
#'@examples \dontrun{FastqCount(sample)}
#'@examples \dontrun{sapply(manyfiles,FastqCount)}
#'
#'@export

FastqCount <- function(file){
  if(length(file)>1){stop("Multiple files handed to counter, use this function for a single file.")}
  #first we use some regex to get the file type as defined by the suffix
  file.type <- gsub(".*[.](.*$)","\\1",file)
  #then we run two expressions using zcat instead of cat if the file is compressed, each time we output a quater of the total line number
  if(file.type=="gz"){
    lines <- system2("echo",args = paste0("$(cat ",file," | zcat | wc -l)"),stdout=TRUE)
    return(as.numeric(lines)/4)
  }
  if(file.type=="fq" | file.type=="fastq"){
    lines <- system2("echo",args = paste0("$(cat ",file," | wc -l)"),stdout=TRUE)
    return(as.numeric(lines)/4)
  }else{message("File type not recognised")}

  }

