# NGSchool2017 materials

Materials prepared by the instructors of the [#NGSchool2017](https://ngschool.eu/2017). 

## Dependencies

Install prerequesities:
- [bioconda](https://bioconda.github.io/)
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

- bioconda & other dependencies
```bash
# install dependencies - AVOID INSTALLING r-base, as it'll mess up with native R installation
conda install bwa htslib spades quast trimmomatic fastqc gmap blast qualimap star busco tophat bowtie2 gawk igv seqtk glimmer exonerate muscle fasttree mcl trimal augustus homer bedtools bbmap gffutils

# non-bioconda
sudo apt install docker docker-engine docker.io varna blast2 macs
sudo groupadd docker
sudo usermod -aG docker $USER
```

- R, Bioconductor and its packages
```bash
# R - need to add R repo first
echo "deb https://www.stats.bris.ac.uk/R/bin/linux/ubuntu $(lsb_release -c | xargs | cut -f2 -d' ')/" | sudo tee -a /etc/apt/sources.list
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo apt-get update

sudo apt install r-base r-base-dev libcurl4-openssl-dev libxml2-dev libcairo2-dev libxt-dev
sudo R
install.packages("plotly"); install.packages("ggplot2")
source("https://bioconductor.org/biocLite.R") # bioconductor
biocLite('BiocInstaller'); biocLite("ATACseqQC"); biocLite("Diffbifnd"); biocLite("affyPLM"); biocLite("arrayMvout"); biocLite("arrayQualityMetrics"); biocLite("gcrma"); biocLite("hgu133acdf"); biocLite("hgu133plus2.db"); biocLite("simpleaffy")
biocLite("RNAprobR"); biocLite("affy"); biocLite("biomaRt"); biocLite("geneplotter"); biocLite("gplots"); biocLite("limma"); biocLite("sva"); biocLite("Rsamtools"); biocLite("BSgenome.Hsapiens.UCSC.hg19")

## scRNA-seq - you may use it through docker image (but it's huuuuuge!)
install.packages("mvoutlier"); install.packages("statmod"); install.packages("pheatmap"); install.packages("ROCR")
source("https://bioconductor.org/biocLite.R")
biocLite('scater'); biocLite('scran'); biocLite("RUVSeq"); biocLite("pcaMethods"); biocLite("SC3")
biocLite("M3Drop"); biocLite("TSCAN"); biocLite("monocle"); biocLite("destiny");
biocLite("edgeR"); biocLite("DESeq2"); biocLite("MAST"); biocLite("MultiAssayExperiment"); biocLite("SummarizedExperiment")
install.packages("devtools");
devtools::install_github("hemberg-lab/scRNA.seq.funcs"); devtools::install_github("JustinaZ/pcaReduce"); devtools::install_github('satijalab/seurat')
devtools::install_github('jw156605/SLICER'); devtools::install_github("hms-dbmi/scde", build_vignettes = FALSE)
```

## Running exercises

```bash
source /ngschool/2017/.bashrc
```