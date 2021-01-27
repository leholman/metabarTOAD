#' MergeReads
#'
#'This function uses usearch to merge paired end Illumina amplicon data.
#'
#'@param maxdiffs Maximum number of allowed mismatches in the alignment between reads. Increase with greater overlap of forward and reverse sequences.
#'@param pctID Minimum percent identity of alignment between forward and reverse read. Decrease with greater overlap of forward and reverse sequences.
#'@param folderwfiles A path specifying folder containing reads to pair. Default is 1.rawreads.
#'@param folderoutput A path specifying output folder to write merged reads into. Default is 2.mergedreads.
#'@param usearchdest Specify the location of the usearch executable if PATH not configured appropriately.
#'
#'
#'@return None
#'
#'
#'@export


MergeReads <- function(maxdiffs=15,pctID=80,folderwfiles="1.rawreads",folderoutput="2.mergedreads",usearchdest="usearch"){
  filedestination <- folderwfiles
  files <- normalizePath(list.files(filedestination,pattern="R1_001.fastq",full.names = TRUE))
  if(length(list.files(filedestination,pattern="R2_001.fastq"))!=length(files)){stop("Number of forward and reverse reads does not match. Please check files.")}
  loopcounter <- 1
  for (file in files){
    message(paste("Merging sample ",loopcounter,"/",length(files)," current file: ",gsub('^.*/(.*_S[0-9]+_L001_R1_001.fastq)','\\1',file),sep=""))
    gsub(filedestination,folderoutput,file)
    output <- gsub('(^.*/.*)_S[0-9]+_L001_R1_001.fastq','\\1.merged.fastq',gsub(filedestination,folderoutput,file))
    mergearg <- paste("-fastq_mergepairs",file,"-relabel @ -fastq_maxdiffs",maxdiffs,"-fastq_pctid",pctID,"-fastqout",output)
    log <- system2(usearchdest,args=mergearg,stdout = TRUE,stderr = TRUE)
    cat(file="log.txt", log , append=T, sep="\n")
    loopcounter <- loopcounter+1
  }

}
