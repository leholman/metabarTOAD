####Running Pipeline

##Here's how it all stiches toghether.

#We start life in an empty folder. Ready to run the analysis

#First we set up the folders required to run the analysis
Folders()

#Great, lets check the folders in finder. We should see 8 folders with different names prefixxed with a number.
list.files()

#We now move our metadata and raw data into the folders. ACTION

#Now we can Unzip our files
Unzip()

#Lets check if our files have been UnZipped
list.files("1.rawreads/", pattern=".fastq")

#Now lets get our new names
index <- read.csv("metadata.csv")

#We now put the desired and current names into the epxression RenameFiles to change the names of our samples
RenameFiles(SeqIDs = as.character(index$RunID),DesiredIDs = as.character(index$RealID))
#We can also rename the samples using an index sheet
RenameFiles(usemetadatafile = TRUE)


#Lets count the number of reads per sample
#gather all the files into one
files <- list.files("1.rawreads",pattern=".fastq",full.names = TRUE)
#now count using FatsqCount
rawreadcount <- sapply(files,FastqCount)

#What is the mean and standard deviation of reads per sample
mean(rawreadcount)
sd(rawreadcount)


#Now we can merge the paired files
MergeReads(usearchdest = "/Users/Luke/Bioinformatics/PATH/usearch")

#Lets count the merged reads
sapply(list.files("2.mergedreads",pattern=".fastq",full.names = TRUE),FastqCount)

#Great now we can stirp the primers off and length truncate






