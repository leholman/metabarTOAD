#'ApplyLulu
#'
#'This function automates the lulu algorithm on a user supplied set of sequences and sample-OTU mapping.
#'
#'
#'@param seqs The location of a FASTA file containing sequences used in the generation of the OTU-sample table.
#'@param table The location of a .csv sample table to be used for lulu OTU pruning. Must contain the same OTU name as the supplied FASTA seqs.
#'@param output An output file location and name for the lulu pruned OTU-sample table .csv file.
#'@param vsearchdest Specify the location of the vsearch executable if PATH not configured appropriately.
#'
#'
#'@return None
#'
#'@importFrom lulu lulu
#'@importFrom utils read.csv write.table read.table
#'
#'@export


ApplyLulu <- function(seqs,table,output,vsearchdest="vsearch",minimum.match=98){
  if(!grepl(".csv",output)){stop("Supplied output is not a .csv file. specify .csv at the end of the file")}
  if(!file.exists("8.LULU")){stop("8.LULU folder not found, please create a folder named 8.LULU in the working directory.")}
  vsearcharg <- paste("--usearch_global ", seqs, " --db ",seqs," --self --id .84 --iddef 1 --userout 8.LULU/",gsub(".*/(.*).csv","\\1",output),".match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10",sep="")
  system2(vsearchdest,vsearcharg)
  rawtable <- read.csv(table)
  OTUselfhits <- read.table(paste("8.LULU/",gsub(".*/(.*).csv","\\1",output),".match_list.txt",sep=""),colClasses = "character")
  OTUselfhits$V3 <- as.numeric(OTUselfhits$V3)
  #extra step here to get rid of all zero rows in the OTU tab
  row_sub = apply(rawtable, 1, function(row) all(row ==0 ))
  rawtable <-rawtable[!row_sub,]
  curatedresults <- lulu(rawtable,OTUselfhits,minimum_match=minimum.match)
  curatedtable <- curatedresults$curated_table
  write.table(curatedtable,output,sep=",")
  ##LULU makes a log file - think of a way to move or get rid of this.
}


