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

4) Kallisto
Kallisto is a program for quantifying transcript abundances through pseudoalignment for
rapid determination of compatibility of reads. Kallisto is also used to build an index from a
reference genome.

