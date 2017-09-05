# NGSchool2017 materials

Materials prepared by the instructors of the [#NGSchool2017](https://ngschool.eu/2017). 

Prerequesities:
- [bioconda](https://bioconda.github.io/)
- [Docker](https://docker.com/)



```bash
# install conda ie. in /ngschool/src/miniconda2
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
chmod +x Miniconda2-latest-Linux-x86_64.sh
./Miniconda2-latest-Linux-x86_64.sh

# configure bioconda channels
(conda config --add channels r)
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```


```bash
# install dependencies
conda install bwa htslib spades quast trimmomatic fastqc gmap blast qualimap star busco tophat bowtie2 gawk igv r r-plotly bioconductor-rnaprobr bioconductor-affy bioconductor-biomart bioconductor-geneplotter r-gplots bioconductor-limma bioconductor-sva seqtk
conda install glimmer exonerate muscle fasttree mcl trimal augustus homer macs2 bioconductor-rsamtools bedtools bioconductor-bsgenome.hsapiens.ucsc.hg19

# non-bioconda
sudo apt install docker docker-engine docker.io varna blast2 (blast-legacy)

# R & bioconductor
sudo apt install R
sudo R
source("https://bioconductor.org/biocLite.R")
biocLite("RNAprobR"); biocLite("affy"); biocLite("biomaRt"); biocLite("geneplotter"); biocLite("gplots")


```


```bash
# import miniconda
export PATH=/ngschool/src/miniconda2/bin:$PATH
```