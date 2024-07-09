(1) ###### Raw reads download from ENA (ENA)
###### open terminal
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR930/008/SRR9304738/SRR9304738.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR930/000/SRR9304740/SRR9304740.fastq.gz

(2)#### trimming (BBduk)
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

 (3)####  kallisto 
###Kallisto is a program for quantifying transcript abundances through pseudoalignment for rapid determination of compatibility of reads. Kallisto is also used to build an index from a reference genome.
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

###Building an indexl
#Next, build an index type:
 ###first download transcript fasta file - https://www.gencodegenes.org/mouse/
kallisto index -i transcripts.idx transcripts.fasta.gz
use cammand - kallisto index -i transcripts.idx gencode.vM35.transcripts.fa (create a transcripts.idx file)

 (4)#########   Quantification using kallisto

##Single end reads
##If your reads are single end only you can run kallisto by specifying the --single flag,


kallisto quant -i transcripts.idx -o output -b 100 --single -l 180 -s 20 reads_1.fastq.gz

use cammand - kallisto quant -i transcripts.idx -o output -b 100 --single -l 180 -s 20 trim38.fq.gz
###paired end reads
kallisto quant -i transcripts.idx -o output -b 100 reads_1.fastq.gz reads_2.fastq.gz

###You can also call kallisto with

kallisto quant -i transcripts.idx -o output -b 100 <(gzcat reads_1.fastq.gz) <(gzcat reads_2.fastq.gz)

(5)###  sleuth

#sleuth is a program for differential analysis of RNA-Seq data. It makes use of quantification uncertainty estimates obtained via kallisto for accurate differential analysis of isoforms or genes, allows testing in the context of experiments with complex designs, and supports interactive exploratory data analysis via sleuth live. The sleuth methods are described in H Pimentel, NL Bray, S Puente, P Melsted and Lior Pachter, Differential analysis of RNA-seq incorporating quantification uncertainty, Nature Methods (2017), advanced access.Scripts reproducing all the results of the paper are available here.
Create a single TSV file that has the TPM abundance estimates for all six samples.

cd $RNA_HOME/expression/kallisto
paste abundance.tsv | cut -f 1,2,5,10,15,20,25,30 > transcript_tpms_all_samples.tsv
ls -1 abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\tlength$_\n"' > header.tsv
cat header.tsv transcript_tpms_all_samples.tsv | grep -v "tpm" > transcript_tpms_all_samples.tsv2
mv transcript_tpms_all_samples.tsv2 transcript_tpms_all_samples.tsv
rm -f header.tsv
paste abundance.tsv | cut -f 1,2,5,10,15,20,25,30 > transcript_tpms_all_samples.tsv




