<img src="/images/Metabartoad.png" align="right" height="60">

# metabarTOAD 

## Introduction
This is a simple pipeline for analysing metabarcoding data. It strings together functions provided in several packages and can be used for the analysis of eDNA metabarcoding data. 

## Requirements
The pipeline has been tested on Ubuntu/Debian and Mac OS 10.14 onwards and (depending on the functions used) requires the following software packages.

* [USEARCH](https://www.drive5.com)
* [VSEARCH](https://github.com/torognes/vsearch) 
* [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html)
* [DADA2](https://benjjneb.github.io/dada2/dada-installation.html)
* [LULU](https://github.com/tobiasgf/lulu)

## Installation
* Install above dependancies 
* Launch R 
```
install.packages("devtools")
library("devtools")
install_github("leholman/metabarTOAD")
library("metabarTOAD")
```
## Pipeline Instructions
You can read an overview of the pipeline using the [UPARSE/UNOISE](http://www.drive5.com/uparse/) algorithm [here](https://github.com/leholman/metabarTOAD/tree/master/vignettes/workflowdada2.pdf), and [DADA2](https://benjjneb.github.io/dada2/index.html) [here](https://github.com/leholman/metabarTOAD/tree/master/vignettes/UPARSEworkflow.pdf).

## Workflow Diagrams
#### UPARSE/UNOISE 
<img src="images/workflow.png" width=600>

#### DADA2 
<img src="images/dada2workflow.png" width=600>
