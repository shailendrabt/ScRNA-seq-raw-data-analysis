# ScRNA-seq-raw-data-analysis


This chapter describes a pipeline for basic bioinformatics analysis of single-cell sequencing data. Starting with raw sequencing data, we will describe how to quality check samples (FastQC), to create an index from a reference genome, to align the sequences to an index, and to quantify transcript abundances. The curated data sets will enable differential expression analysis, population analysis, and pathway analysis.

1)Download Raw reads
first download raw data from ENA OR SRA or GEO in fastq.gz format.
###### open terminal
wget  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR930/008/SRR9304738/SRR9304738.fastq.gz

wget  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR930/000/SRR9304740/SRR9304740.fastq.gz

2) FastQC
When first obtaining raw sequencing data, check the quality of the sequencing experiment.
Installation of this program does not require the use of command line. This program requires
Java 7 or higher to run. Website to download FASTQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.

3) BBDuk
  BBDuk is a trimming tool that is part of the BBTools toolset. BBTools has a variety of other
bioinformatics tools; however, here only BBDuk is used. Download it from sourceforge.net/
projects/bbmap/. Move the downloaded file to a working directory of choice.
To install, go to the terminal or shell and change to the working directory
and then input the following command:

$ wget https://github.com/BioInfoTools/BBMap/releases/download/v35.85/BBMap_35.85.tar.gz

$ tar --xzvf BBMap_(version).tar.gz ( version- BBMap_35.85.tar.gz )
#### trimming (BBduk)
##Adapter trimming:
bbduk.sh -Xmx1g in=SRR9304738.fastq.gz out=trim38.fq.gz qtrim=rl trimq=10
bbduk.sh -Xmx1g in1=read1.fq in2=read2.fq out1=clean1.fq out2=clean2.fq ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo

###USE CAMMAND 
bbduk.sh -Xmx1g in=SRR9304740.fastq.gz out=trim40.fq.gz qtrim=rl trimq=10
 bbduk.sh -Xmx1g in=SRR9304738.fastq.gz out=trim38.fq.gz qtrim=rl trimq=10
 
###Quality trimming:
bbduk.sh -Xmx1g in=reads.fq out=clean.fq qtrim=rl trimq=10
bbduk.sh -Xmx1g in=read1.fq in=read2.fq out1=clean1.fq out2=clean2.fq qtrim=rl trimq=10

####Contaminant filtering:
bbduk.sh -Xmx1g in=reads.fq out=unmatched.fq outm=matched.fq ref=phix.fa k=31 hdist=1 stats=stats.txt
or
bbduk.sh -Xmx1g in1=r1.fq in2=r2.fq out1=unmatched1.fq out2=unmatched2.fq outm1=matched1.fq outm2=matched2.fq ref=phix.fa k=31 hdist=1 stats=stats.txt

###Kmer masking:
bbduk.sh -Xmx1g in=ecoli.fa out=ecoli_masked.fq ref=salmonella.fa k=25 ktrim=N rskip=0

4) Kallisto
Kallisto is a program for quantifying transcript abundances through pseudoalignment for
rapid determination of compatibility of reads. Kallisto is also used to build an index from a
reference genome.
##open terminal & paste cammands
##Download and installation
git clone https://github.com/pachterlab/kallisto
cd kallisto
mkdir build
cd build
cmake ..
make
make install
kallisto 
kallisto index
kallisto quant
or 
To install kallisto, use the following commands in shell:
$ ruby --e “$(curl --fsSL https://raw.githubusercontent.com/Homebrew/install/
master/install)”
$ brew tap homebrew/science
Author Manuscript
$ brew install kallisto

5) Sleuth
Sleuth is an R package designed for usage downstream of Kallisto.
To install Sleuth use the following packages in R:
> source(“http://bioconductor.org/biocLite.R”)
> biocLite(“rhdf5”)
> install.packages(“devtools”)
> devtools::install_github(“pachterlab/sleuth”)


6) Singular
Singular is an analysis toolset offered gratis through Fluidigm. Its ease of use through a
graphical interface makes it easy for beginners in bioinformatics to perform basic analysis.
However, this R package does not work with Mac computers. Singular analysis toolset
software and practice sets can be downloaded from www.fluidigm.com/software. To install
the package, in R, go to the packages tab and click install packages from local files. Then,
select the downloaded zipped file.
Once installed, in R, type:
> Library(fluidigmSC)
> firstRun()

7) SCDE
SCDE is an R package used in the statistical analysis of single-cell RNA-seq data. Use this
to observe differential expression across samples.
To install SCDE use the following packages in R:
> source(“https://bioconductor.org/biocLite.R”)
Author Manuscript
> biocLite(“scde”)

8) DAVID
To use DAVID, proceed to the following URL: david.ncifcrf.gov.
