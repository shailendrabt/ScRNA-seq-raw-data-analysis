3) ####    kallisto 
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
###Building an index\
#Next, build an index type:
 ###first download transcript fasta file - https://www.gencodegenes.org/mouse/
kallisto index -i transcripts.idx transcripts.fasta.gz
use cammand - kallisto index -i transcripts.idx gencode.vM35.transcripts.fa (create a transcripts.idx file)

#########Quantification

##Single end reads
##If your reads are single end only you can run kallisto by specifying the --single flag,


kallisto quant -i transcripts.idx -o output -b 100 --single -l 180 -s 20 reads_1.fastq.gz
use cammands- kallisto quant -i transcripts.idx -o output -b 100 --single -l 180 -s 20 trim38.fq.gz
            - kallisto quant -i transcripts.idx -o output -b 100 --single -l 180 -s 20 trim40.fq.gz

###paired end reads
kallisto quant -i transcripts.idx -o output -b 100 reads_1.fastq.gz reads_2.fastq.gz
