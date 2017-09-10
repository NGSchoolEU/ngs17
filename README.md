# NGSchool2017 materials

Materials prepared by the instructors of the [#NGSchool2017](https://ngschool.eu/2017). 

**Table of Contents**  
   * [NGSchool2017 materials](#ngschool2017-materials)
      * [Dependencies](#dependencies)
      * [Running exercises](#running-exercises)
         * [Working in your own laptop](#working-in-your-own-laptop)
         * [Working in remote NGSchool server](#working-in-remote-ngschool-server)
         * [Working in VirtualBox](#working-in-virtualbox)
      * [Cloning the repository](#cloning-the-repository)
         * [Materials not included in github repo](#materials-not-included-in-github-repo)
            * [Introduction](#introduction)
            * [De novo assembly](#de-novo-assembly)    
            * [Hi-C](#hi-c)    
            * [Microbial genomics](#microbial-genomics)


## Dependencies
In order to run workshop examples in your own laptop, you'll need to install all below prerequesities.  
**Note, the installation instructions are meant for Ubuntu 16.04. 
Everything should be done in below order, it may take 2-3 hours and around 15-20GB of hard-drive space.**

### [bioconda](https://bioconda.github.io/) & [docker](https://docker.com)
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

# install dependencies - AVOID INSTALLING r-base, as it'll mess up with native R installation
conda install bwa htslib samtools spades quast trimmomatic fastqc gmap blat blast qualimap star busco tophat bowtie2 gawk igv seqtk glimmer exonerate muscle fasttree mcl trimal augustus homer bedtools bbmap gffutils

# non-bioconda
sudo apt install docker.io varna blast2 macs
sudo groupadd docker
sudo usermod -aG docker $USER
```

### R, Bioconductor and other R packages
```bash
# R - need to add R repo first
echo "deb https://www.stats.bris.ac.uk/R/bin/linux/ubuntu $(lsb_release -c | xargs | cut -f2 -d' ')/" | sudo tee -a /etc/apt/sources.list
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo apt-get update && sudo apt upgrade && sudo apt install libcurl4-openssl-dev libxml2-dev libcairo2-dev libxt-dev libssl-dev
sudo apt install r-base r-base-dev

# install R packages for all users
sudo R
install.packages("plotly"); install.packages("ggplot2")
source("https://bioconductor.org/biocLite.R") # bioconductor
biocLite('BiocInstaller'); biocLite("ATACseqQC"); biocLite("Diffbifnd"); biocLite("affyPLM"); biocLite("arrayMvout"); biocLite("arrayQualityMetrics"); biocLite("gcrma"); biocLite("hgu133acdf"); biocLite("hgu133a.db"); biocLite("hgu133plus2.db"); biocLite("simpleaffy")
biocLite("RNAprobR"); biocLite("affy"); biocLite("biomaRt"); biocLite("geneplotter"); biocLite("gplots"); biocLite("limma"); biocLite("sva"); biocLite("Rsamtools"); biocLite("ChIPseeker"); 
biocLite("BSgenome.Hsapiens.UCSC.hg19") # large

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

### 8/09/2017 UPDATE
```bash
sudo R # tgambin & kkedzierska
install.packages("data.table");
source("https://bioconductor.org/biocLite.R")
biocLite('parallel'); biocLite('RCurl'); biocLite('gdata'); biocLite('Hmisc'); biocLite('matrixStats'); biocLite('DNAcopy'); biocLite('GenomicRanges'); biocLite('Rsubread'); biocLite('WES.1KG.WUGSC'); biocLite('CODEX'); biocLite("ChIPseeker");
```

### Manual installation
  - [rnaQUAST](http://cab.spbu.ru/software/rnaquast/)

## Running exercises

### Working in your own laptop
Copy workshop materials locally ie. `rsync -av /media/$USER/USB_MOUNT_DIR ~/ngschool/2017`,
enter NGSchoool directory `cd ~/ngschool/2017` and you are ready to work.
Make sure, you have installed [all prerequesities](#dependencies) before! 

### Working in remote NGSchool server
Login to the server with your credentials,
sync workshop materials to your home directory `rsync -av --exclude '*.git/' /ngschool/2017 ~/ngschool`, 
enter your personal NGSchoool directory `cd ~/ngschool/2017` and you are ready to work.

Make sure to import local variable in each new window
```bash
# NOTE: you may need to change `/ngschool/2017` directory
# if you cloned the repository to another location
source /ngschool/2017/.bashrc
```

### Working in VirtualBox
First, [get VM image](http://zdglab.iimcb.gov.pl/cluster/ngschool/2017/VM/Ubuntu16.04.3.vdi)
and [create VM in VirtualBox using this image](http://linuxbsdos.com/2015/11/13/how-to-import-a-virtual-machine-image-into-virtualbox/). 
Then run VM (u: ngschool p: ngschool), enter NGSchoool directory `cd /ngschool/2017` and you are ready to work. 

## Cloning the repository
**This has to be done only if you wish to explore materials before the school. Otherwise, ignore below.**  
To clone repo, use `git clone --recursive https://github.com/NGSchoolEU/2017.git`  
Below, we're providing links to data not included in this repository. 

### Materials not included in github repo
You can get below using `wget -nc -r -np HTTP`

#### Introduction
All exercises are in: http://compbio.fmph.uniba.sk/temp/ngschool2017/
- http://compbio.fmph.uniba.sk/temp/ngschool2017/everything.zip

#### De novo assembly
- http://spades.bioinf.spbau.ru/~school/genomics/
- http://spades.bioinf.spbau.ru/~school/transcriptomics/

#### Hi-C
Get ONE of below files: 
- http://makarich.fbb.msu.ru/agalicina/2017/NGS/hic_workshop_2017_withimage.tar.gz (2GB, docker image compiled)
- http://makarich.fbb.msu.ru/agalicina/2017/NGS/hic_workshop_2017_light.tar.gz (0.7GB, without docker image compiled)

#### Microbial genomics
- http://bugs.sgul.ac.uk/ngschool-2017/
