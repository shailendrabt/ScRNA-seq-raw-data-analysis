# ScRNA-seq-raw-data-analysis


This chapter describes a pipeline for basic bioinformatics analysis of single-cell sequencing data. Starting with raw sequencing data, we will describe how to quality check samples (FastQC), to create an index from a reference genome, to align the sequences to an index, and to quantify transcript abundances. The curated data sets will enable differential expression analysis, population analysis, and pathway analysis.

A)Download Raw reads     #####################################################################################

first download raw data from ENA OR SRA or GEO in fastq.gz format.
###### open terminal
wget  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR930/008/SRR9304738/SRR9304738.fastq.gz

wget  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR930/000/SRR9304740/SRR9304740.fastq.gz

#QC and Trimming of Raw Sequences

Load fastq files into FASTQC.
2.Observe per base sequence quality. If the curve dips below a Phred score of 28,
the sequence will require trimming (see Note 1). If more than half of the
sequence is below 28, file may not be adequate for downstream analysis and
should be removed. Any adapter sequences present in the overrepresented
sequences will also need to be removed.
3.Perform the following command in terminal or shell. For the inputs, direct the
pathway to the location of files. For pairedend sequencing data, provide two
inputs, one for each sample. The outputs will be designated with a new directory.
The qtrim parameter designates the direction of the trim (either left or right). In
the following example, right is selected with “r.” The trimq parameter designates
the Phred score threshold where trimming should occur. Here we set the
threshold to 28, which will trim the bases to the right of where the Phred score
dips below 28.

B) FastQC             #########################################################################################
When first obtaining raw sequencing data, check the quality of the sequencing experiment.
Installation of this program does not require the use of command line. This program requires
Java 7 or higher to run. Website to download FASTQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.

C) BBDuk   ######################################################################################################
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

D) Kallisto  ###########################################################################################
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
###Index Building
Perform the following command in the terminal or shell. The generated index file
will be output into the current directory. Change directory prior to this command.
##Alignment and Quantification
1.Perform the following command after building the index. Here we set the
number of bootstraps to 100. This can be changed if desired. The input files are
the trimmed fastq inputted as pairs.
To perform pseudoalignment:
$ kallisto quant --i /path/to/nameOfIndex -b 100 -o nameOfDirectoryOutput /
path/to/input1_1.fastq / path/to/input1_2.fastq
2.To expedite the process of quantifying the samples, write a script to automate the
process. Create a plain text file that lists the name of the fastq files. Then, create
a plain text file with the following commands.
3.Alternatively, use the following command to perform the same function as
“quant” but without bootstrap. First make a plain text file called “batch.txt,”
which includes columns for #id, file1, and file2 names.

E) Sleuth(Filtering Transcript Abundances and Annotation with Sleuth)            #######################################
Sleuth is an R package designed for usage downstream of Kallisto.
1.To begin, run the following commands to open necessary R packages.

To install Sleuth use the following packages in R:
> source(“http://bioconductor.org/biocLite.R”)
> biocLite(“rhdf5”)
> install.packages(“devtools”)
> devtools::install_github(“pachterlab/sleuth”)
2. The following commands will designate the pathway to the location of the
kallisto files. Use .h5 files as they contain boot-straps. First designate the base
directory.
3.Command to designate the directory containing abundance.h5 files.
4.Create a .txt file that lists the name of the samples and their condition. First
column is labeled “Sample_id,” and the second column is “Condition.” This will
be used for the following command.
5. The following commands will now designate the file location to the sample of
interest.
6. The following command sets the path for the s2c variable to the kal_dirs directory.
Go to the summaries tab and to kallisto table to view the transcript abundances.
Save this table as it will be needed to generate a data matrix. Processed data can
also be viewed for alignment QC metrics.
> 
F) Singular          ####################################################################
Singular is an analysis toolset offered gratis through Fluidigm. Its ease of use through a
graphical interface makes it easy for beginners in bioinformatics to perform basic analysis.
However, this R package does not work with Mac computers. Singular analysis toolset
software and practice sets can be downloaded from www.fluidigm.com/software. To install
the package, in R, go to the packages tab and click install packages from local files. Then,
select the downloaded zipped file.
Once installed, in R, type:
> Library(fluidigmSC)
> firstRun()

G) SCDE       ##########################################################################
SCDE is an R package used in the statistical analysis of single-cell RNA-seq data. Use this
to observe differential expression across samples.
To install SCDE use the following packages in R:
> source(“https://bioconductor.org/biocLite.R”)
Author Manuscript
> biocLite(“scde”)

H) DAVID        ###########################################################################
To use DAVID, proceed to the following URL: david.ncifcrf.gov.
