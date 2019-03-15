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
  return(as.numeric(system2("awk",args = paste("'{s++}END{print s/4}'",file),stdout=TRUE)))
  }

