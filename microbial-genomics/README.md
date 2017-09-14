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
* Build an index of the reference genome

```bash
tmap index -f genomes/NC_000962.fna
```

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

## Site calling

We now need to process the alignment file and identify which base is present at every reference genome position, thus identifying mutations or wild type at sites of interest. We use the “samtools mpileup” program to process the alignment file and call each genome site where the reads have mapped [type all on one line]:

```bash
samtools-mpileup genomes/NC_000962.fna Case1a.sort.bam > Case1a.all.vcf.gz
```

This will take a little time, so we will move on to the next stage and come back to this when it is finished

## Visualise alignment file using Artemis

Artemis is a genome visualisation tool developed at the Sanger Institute. It is available for download from [http://www.sanger.ac.uk/science/tools/artemis]

* Artemis can load a genome file in GenBank format and overly the bam alignment file from the command line [Note, ignore the warnings by clicking No in the popup]:

```bash
java -Dbam=Case1a.sort.bam -jar bin/artemis.jar genomes/NC_000962.gbk
```

* Take some time to navigate around the interface
* Right click on the alignment panel and click Show->SNP marks, this will show SNP positions in read
* Navigate to genome base position 7570, by clicking Goto->Navigator and entering the position in "Goto Base:", then click Goto.
* What effect does the mutation have in Case1a?

## Site calling (continued)

When the "samtools mpileup" command  has finished we can continue.

* build a tabix index, (tabix needs this below)

```bash
tabix -p vcf Case1a.all.vcf.gz
```

We can then filter the resulting VCF (Variant Calling Format) file and examine those sites known to be associated with drug resistance phenotypes.

Identify which amino acids have been changed below, a codon translation table is available here [https://www2.le.ac.uk/projects/vgec/diagrams/34%20codon%20table.jpg].

1. Rifampicin
* RpoB
* amino acid 430, L -> P

```bash
tabix vcf/Case1a.all.vcf.gz NC_000962.3:761094-761096
```

* amino acid 450, S -> L

```bash
tabix vcf/Case1a.all.vcf.gz NC_000962.3:761154-761156
```

2. Isoniazid
* KatG
* amino acid 315, S -> T [Note that katG is coded on the reverse strand]

```bash
tabix vcf/Case1a.all.vcf.gz NC_000962.3:2155167-2155169
```

3. Ethambutol
* EmbB
* amino acid  306, M -> I / V

```bash
tabix vcf/Case1a.all.vcf.gz NC_000962.3:4247429-4247431
```

* amino acid 406, G -> A

```bash
tabix vcf/Case1a.all.vcf.gz NC_000962.3:4247729-4247731
```

4. Fluoroquinolones
* GyrA
* amino acid 90, A -> V

```bash
tabix vcf/Case1a.all.vcf.gz NC_000962.3:7569-7571
```

* amino acid 94, D -> G / E

```bash
tabix vcf/Case1a.all.vcf.gz NC_000962.3:7581-7583
```

Try the analysis process for other Cases by changing the filename Case1a.all.vcf.gz above, compare to the attached published results for this set of isolates. To get a full list of patient sample files type:

```bash
ls vcf/
```

## Phylogenetic analysis

Phylogenetic reconstruction is used to assess the level of similarity between the genomes of each isolate, and places them into the context of an evolutionary tree.
As above VCF files are generated for all isolates to be examined. These VCF’s are then filtered at each site according to several criteria and any site failing is removed from the analysis. The criteria are:

* Quality (Q) score, > 30
* Depth of coverage, at least 4 reads required
* More than 75% reads support site call
* Site present in all sequences (ignores insertion sequences etc)

The remaining sites are collected into a matrix
To perform this site filtering and matrix building (this will take a few minutes to complete)

```bash
snp-filter
```

Phylogenetic reconstruction is performed by Maximum Likelihood estimation, using RAxML

```bash
estimate-tree
```

Visualise tree file using FigTree
* Download FigTree from http://tree.bio.ed.ac.uk/software/figtree/
* Double Click to decompress the downloaded FigTree file
* Download tree file to Desktop (RAxML_bipartitions.NC_000962.SGUL.b1)
* From FigTree click Open and find your tree file (RAxML_bipartitions.NC_000962.SGUL.b1)
* Explore the tree with the various FigTree options e.g. root the tree at the midpoint


## Automated tools

Above we have seen the “manual” approach to resistance prediction in TB. For WGS to be useful tool in clinical practice, these analysis pipelines need to be standardised and packaged into automated tools Several such tools are available, try uploading the fastq files in the data directory (you will need to download the fastq.gz file to your PC first using the CyberDuck as above) into these tools:

1. TB profiler
* Developed by Prof Taane Clark’s group at LSHTM.
* Open a web browser and load http://tbdr.lshtm.ac.uk/
* Upload Case1a.fastq.gz using the “Gzipped fastq file” “Choose file” button.
* Click Submit
* Compare results with those determined in the workshop

2. Mykrobe predictor
* Developed by Dr Zamin Iqbal at University of Oxford.
* Open a web browser and load http://www.mykrobe.com/products/predictor/
* Select the TB tab and click Download
* Agree to the licence and Download the version for your OS.
* Run your version of the software, this depends on which OS you are using
* Drag the Case1a.fastq.gz file to the Predictor window, or click Browse to find.
* Explore the reports and compare the results as above


