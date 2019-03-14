#' UnZip Files
#'
#'This function unzips gzip files. Illumina currently zips all raw sequencing
#'output in gzip. The fucntion looks for gzip files in the 1.rawreads folder,
#'alternatively another folder can be sepcified.
#'It relies on pigz so rock on over to 'https://zlib.net/pigz/' and
#'install a copy.
#'
#'@param folderwfiles Specify a folder other than '1.rawreads' to look for Gzip files.
#'@param unpigzlocation If the 'unpigz' executable in not in your path specify the location here.
#'
#'@return None
#'
#'@examples Unzip()
#'
#'@export


Unzip <- function(folderwfiles="1.rawreads",unpigzlocation=""){
  files <- list.files(folderwfiles, pattern=".gz", full.names = TRUE)
  if(length(files)==0){stop("No gz files found to unzip, check folder.")}
  message(paste("Unzipping",length(files),"files."))
  for (i in 1:length(files)){
    system2(paste(unpigzlocation,"unpigz",sep=""),args=files[i])
  }
}


