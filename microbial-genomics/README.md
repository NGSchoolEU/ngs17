# Workshop - Drug resistance phenotype prediction using WGS

In this workshop we will use sequence data generated on an IonTorrent Personal Genome Machine to identify any known mutations associated with drug resistance, therefore predicting the drug resistance phenotype. We will use standard tools on linux to perform the analyses.

## Setup files

* If you don't already have the data then rsync it using

```bash
rsync -avz --exclude="*.git/" USERNAME@192.168.1.111:/ngschool/2017/microbial-genomics ~/ngschool/
```

Then change into the working directory

```bash
cd ~/ngschool/microbial-genomics
```
* add the local bin directory to your PATH

```bash
export PATH=./bin:$PATH
```

## Data files

* Data is available for 16 isolates from 6 suspected XDR-TB cases. The study is described in “Clinical application of whole-genome sequencing to inform treatment for multidrug-resistant tuberculosis cases. J Clin Microbiol. 2015 May;53(5):1473-83. doi: 10.1128/JCM.02993-14”. [https://www.ncbi.nlm.nih.gov/pubmed/25673793]
* Due to file size constraints only two isolate sequences are included initially, but the remaining sequence files are available from the local workshop server (see below).
* We will use WGS data to predict drug resistance phenotypes and therefore inform treatment options. 
* Explore the data files here:

```bash
ls data/
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
tmap mapall -f genomes/NC_000962.fna -r data/Case1a.bam -n 1 -v -Y -u -o 1 stage1 map4 > Case1a.bam
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

We now need to process the alignment file and identify which base is present at every reference genome position, thus identifying mutations or wild type at sites of interest. We use the “samtools mpileup” program to process the alignment file and call each genome site where the reads have mapped [enter all on one line]:

```bash
samtools mpileup -ugf genomes/NC_000962.fna Case1a.sort.bam | bcftools view -cg - | bgzip > Case1a.all.vcf.gz
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
1. Which gene is this mutation found in?
2. What effect does the mutation have in Case1a? (Check if this is a known TB resistance mutation by checking the TBDreamDB database [https://tbdreamdb.ki.se/Data/MutationDetail.aspx?AreaId=FLQ&GeneID=Rv0006&OnlyHC=true])

## Site calling (continued)

When the "samtools mpileup" command  has finished we can continue.

* build a tabix index, (tabix needs this below)

```bash
tabix -p vcf Case1a.all.vcf.gz
```

We can then filter the resulting VCF (Variant Calling Format) file and examine those sites known to be associated with drug resistance phenotypes.

Identify which amino acids have been changed below, a codon translation table is available here [https://www2.le.ac.uk/projects/vgec/diagrams/34%20codon%20table.jpg].

In these vcf formated files the second column shows the genome position, the fourth column is the reference base while the fifth column shows if there is a SNP at that position, if there is no SNP then a '.' is shown. These commands show the three bases of the relevant codon.

#### Rifampicin
Resistance to Rifampcicin is driven by mutations in the RNA polymerase B gene, *rpoB*. We will look at two common sites here:

* RpoB - amino acid 430, L -> P

```bash
tabix vcf/Case1a.all.vcf.gz NC_000962.3:761094-761096
```

* RpoB - amino acid 450, S -> L

```bash
tabix vcf/Case1a.all.vcf.gz NC_000962.3:761154-761156
```

Is Case1a resistant to Rifampicin?

#### Isoniazid
Common mutations to Isoniazid are found in the *inhA* and *katG* genes. We will look at a common *katG* mutation here:

* KatG - amino acid 315, S -> T [Note that katG is coded on the reverse strand]

```bash
tabix vcf/Case1a.all.vcf.gz NC_000962.3:2155167-2155169
```
Is Case1a resistant to Isoniazid?

#### Ethambutol
Mutations in the *embB* gene are commonly associated with reistsnace to Ethambutol.

* EmbB - amino acid  306, M -> I / V

```bash
tabix vcf/Case1a.all.vcf.gz NC_000962.3:4247429-4247431
```

* EmbB - amino acid 406, G -> A

```bash
tabix vcf/Case1a.all.vcf.gz NC_000962.3:4247729-4247731
```

Is Case1a resistant to Ethambutol?

#### Fluoroquinolones
Resistance to fluoroquinolones is assocated with mutations in the A subunit of the DNA gyrase gene, *gyrA*, we will look at two common sites here, these are located within the Quinolone Resistance determining Region (QRDR):

* GyrA - amino acid 90, A -> V (causes low level resistance to moxifloxacin)

```bash
tabix vcf/Case1a.all.vcf.gz NC_000962.3:7569-7571
```

* GyrA - amino acid 94, D -> G / E (causes high level resistance to moxifloxacin)

```bash
tabix vcf/Case1a.all.vcf.gz NC_000962.3:7581-7583
```

Is Case1a resistant to Fluoroquinolones?

Compare your results to the published results for this set of isolates found here [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4400773/table/T2/]

Try the analysis process for other Cases by changing the filename Case1a.all.vcf.gz above. To get a full list of patient samples list the vcf files available:

```bash
ls vcf/
```

Or you can rerun the analysis using the other isolate raw data files. Case3a is also provided and the remaining data files can be obtained using:

```bash
rsync -avz --exclude="*.git/" USERNAME@192.168.1.111:/ngschool/users/awitney/data/* ~/ngschool/mcirobial-genomics/data/
```

## Phylogenetic analysis

Phylogenetic reconstruction is used to assess the level of similarity between the genomes of each isolate, and places them into the context of an evolutionary tree.
As above VCF files are generated for all isolates to be examined. These VCF’s are then filtered at each site according to several criteria and any site failing is removed from the analysis. The criteria are:

* Quality (Q) score, > 30
* Depth of coverage, at least 4 reads required
* More than 75% reads support site call
* Site present in all sequences (ignores insertion sequences etc)

The remaining sites are collected into a matrix
To perform this site filtering and matrix building (this will take a few minutes to complete), use the following command [note, this is also contained in the snp-filter script] 

```bash
perl scripts/snp_caller.pl --chrom NC_000962.3 --qual 30 --dp 2 --dp4 75 --dpmax 5000 --af 0 --mq 30 --noindels --noheader \
  --vcf Case1a,vcf/Case1a.all.vcf.gz \
  --vcf Case2a,vcf/Case2a.all.vcf.gz \
  --vcf Case3a,vcf/Case3a.all.vcf.gz \
  --vcf Case4a,vcf/Case4a.all.vcf.gz \
  --vcf Case5a,vcf/Case5a.all.vcf.gz \
  --vcf Case6a,vcf/Case6a.all.vcf.gz \
  --vcf Case7a,vcf/Case7a.all.vcf.gz \
  --vcf Case8a,vcf/Case8a.all.vcf.gz \
  --dir ./ --phylip NC_000962.b1.infile --verbose 0 --refilter -b 1 --cpus 1
```

Phylogenetic reconstruction is performed by Maximum Likelihood estimation, using RAxML

```bash
raxmlHPC-PTHREADS-SSE3 -T 1 -f a -s NC_000962.b1.infile -x 12345 -p 1234 -# 100 -m GTRGAMMA -n NC_000962.b1 -o NC_000962.3
```

Visualise tree file using FigTree
* Download FigTree from http://tree.bio.ed.ac.uk/software/figtree/
* Double Click to decompress the downloaded FigTree file
* Download tree file to Desktop (RAxML_bipartitions.NC_000962.SGUL.b1)
* From FigTree click Open and find your tree file (RAxML_bipartitions.NC_000962.SGUL.b1)
* Explore the tree with the various FigTree options e.g. root the tree at the midpoint


## Automated tools

Above we have seen the “manual” approach to resistance prediction in TB. For WGS to be useful tool in clinical practice, these analysis pipelines need to be standardised and packaged into automated tools Several such tools are available, try uploading the fastq files in the data directory (you will need to download the fastq.gz file to your environment/PC first) into these tools:

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


