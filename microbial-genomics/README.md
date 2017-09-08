# Workshop - Drug resistance phenotype prediction using WGS

In this workshop we will use sequence data generated on an IonTorrent Personal Genome Machine to identify any known mutations associated with drug resistance, therefore predicting the drug resistance phenotype. Data is provided for 16 isolates from 6 suspected XDR-TB cases. We will use standard tools on linux to perform the analyses.

## Setup and server login

* First change to the workshop directory

```bash
cd microbial-genomics
```

## Data files

* Data is provided for 16 isolates from 6 suspected XDR-TB cases. The study is described in “Clinical application of whole-genome sequencing to inform treatment for multidrug-resistant tuberculosis cases. J Clin Microbiol. 2015 May;53(5):1473-83. doi: 10.1128/JCM.02993-14”
* We will use WGS data to predict drug resistance phenotypes and therefore inform treatment options. 
* Explore the data files here:

```bash
ll data/
```

* The Ion Torrent platform generates raw data in unaligned bam files rather than fastq. Take a look inside a raw data file, you will notice that it looks exactly the same as an aligned bam, but without the reference and coordinates

```bash
samtools view –h data/Case1a.bam | less
```

## Alignment

We must first align the reads generated in the sequencing run to a standard reference genome. In this case we will use the H37Rv genome available from the RefSeq database at NCBI (accession number NC_000962.3). This process includes the following steps:
* Map the raw reads to the H37Rv reference using TMAP (Note that TMAP is IonTorrent specific, if we were using MiSeq data we would probably use a different aligner e.g. BWA) [Note, type all on one line]

```bash
tmap-mapall genomes/NC_000962.fna data/Case1a.bam > Case1a.bam
```

* Sort the alignment by genome position

```bash
samtools sort Case1a.bam Case1a.sort
```

* Build index

```bash
samtools index Case1a.sort.bam
```

* Assess alignment file for number of reads mapped

```bash
samtools flagstat Case1a.sort.bam
```

* Take a look inside the alignment file [press q to exit] Note, type all on one line] 

```bash
samtools view Case1a.sort.bam | less
```


