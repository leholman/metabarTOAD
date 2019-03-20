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
#We can either provide the primer and length truncation data if we have a single set of primers
PrimerStrip(PrimerF = "NNNNNNGGWACWGGWTGAACWGTWTAYCCYCC",
            PrimerR = "TAIACYTCIGGRTGICCRAARAAYCA",
            MinLen = 303,
            MaxLen = 323,
            cutadaptdest = "/Users/Luke/Bioinformatics/PATH/cutadapt",
            ncores=8)
#or provide data for multiple primer sets and process a bunch toghether
PrimerStrip(UsePrimerFile = TRUE,cutadaptdest = "/Users/Luke/Bioinformatics/PATH/cutadapt",ncores=8)


##Now we need to pool our sequences, get rid of those with errors and fiter out singletons.
#this is done using a single function as below if we used a single primer set for this analysis
PoolNFilterReads(vsearchdest="/Users/Luke/Bioinformatics/PATH/vsearch")

#alternativly if we have a bunch of primers we can run the function like this
PoolNFilterReads(UsePrimerFile=TRUE,vsearchdest="/Users/Luke/Bioinformatics/PATH/vsearch")

#Now if we want to generate 97% clusters and map our reads we can use a single function as below for a single primer set
OTUCluster(usearchdest = "/Users/Luke/Bioinformatics/PATH/usearch")

#and like this if we have used the primer file.
OTUCluster(UsePrimerFile = TRUE,usearchdest = "/Users/Luke/Bioinformatics/PATH/usearch")


##If we prefer actual sequence variants we can generate these in the same way with the below function for a single primer set
DenoiseUNOISE3(usearchdest = "/Users/Luke/Bioinformatics/PATH/usearch")

#and like this if we have used the primer file.
DenoiseUNOISE3(UsePrimerFile = TRUE,usearchdest = "/Users/Luke/Bioinformatics/PATH/usearch")

#If we want to use LULU to curate our data we can do so as below.
ApplyLulu(seqs="5.OTUs/AllSamples.unoise3.OTUs.fasta",
          table="6.mappings/OTUtabs/AllSamplesraw.unoise3.csv",
          output="8.LULU/curateduoise3table.csv",
          vsearchdest="/Users/Luke/Bioinformatics/PATH/vsearch")



