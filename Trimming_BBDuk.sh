BBDuk Guide
“Duk” stands for Decontamination Using Kmers. BBDuk was developed to combine most common data-quality-related trimming, filtering, and masking operations into a single high-performance tool. It is capable of quality-trimming and filtering, adapter-trimming, contaminant-filtering via kmer matching, sequence masking, GC-filtering, length filtering, entropy-filtering, format conversion, histogram generation, subsampling, quality-score recalibration, kmer cardinality estimation, and various other operations in a single pass. Specifically, any combination of operations is possible in a single pass, with the exception of kmer-based operations (kmer trimming, kmer masking, or kmer filtering); at most 1 kmer-based operation can be done in a single pass. BBDuk2 allows multiple kmer-based operations in a single pass, and is otherwise equivalent to BBDuk.
2) #### trimming (BBduk)
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
