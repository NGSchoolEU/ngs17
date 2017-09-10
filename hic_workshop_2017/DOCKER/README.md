### Hi-C data analysis

This is an annotation for the workshop on Hi-C data analysis, NGS'17, Warsaw, Poland. 
Author: Aleksandra Galitsyna agalicina@gmail.com

We will be working in Jupyter notebook environment supported by Docker container. 
For Hi-C data analysis we will use widespread software hiclib, bowtie2 and lavaburst.

Initial data was acquired from public resources (GEO, UCSC) and adapted for laptop computations (e.g. only reads for chromosome 1 were retained). List of reference papers:

- "Single-nucleus Hi-C reveals unique chromatin reorganization at oocyte-to-zygote transition" Nature 544, 110–114 Ilya M. Flyamer,	Johanna Gassler,	Maxim Imakaev,	Hugo B. Brandão,	Sergey V. Ulianov, Nezar Abdennur,	Sergey V. Razin,	Leonid A. Mirny	& Kikuë Tachibana-Konwalski

- "A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping" Cell 159:7, 1665–1680, Suhas S.P. Rao, Miriam H. Huntley, Neva C. Durand, Elena K. Stamenova, Ivan D. Bochkov, James T. Robinson, Adrian L. Sanborn, Ido Machol, Arina D. Omer, Eric S. Lander, Erez Lieberman Aiden

#### Files annotation

Workshop folder for Hi-C analysis consists of 4 main parts:

1. __DATA__ -- folder for data storage. Initially, there are two subdirectories: FASTQ with reads and ANNOT with some pre-computed additional files (such as Hi-C datasets from other papers etc.). During workshop there will be also SAM and HDF5 folders. 

2. __GENOMES__ -- folder with *Homo sapiens* reference genome (HG19_FASTA) and bowtie 2 index (HG19_IND)

3. __WD__ -- folder with Jupyter notebooks that could be run from Docker container. 

4. __DOCKER__ -- folder with this annotation, Dockerfile and compiled Docker image (no image in the light version of this archive!).

#### Docker installation on Ubuntu

[Most useful guide so far](https://docs.docker.com/engine/installation/linux/docker-ce/ubuntu/), in brief (online only):

sudo apt-get update
sudo apt-get install docker-ce
sudo docker run hello-world

#### Docker run from existing image

sudo docker load -i hiclib-notebook.tar

pwd .

sudo docker run -it -v <path_to_workshop_dir>/:/home/jovyan/ -p 8888:8888 hiclib-notebook

#### Docker build for advanced users

If you would like to use Dockerfile or its components, it can be found in the current folder (DOCKER). Compilation of new docker image is dependent on the internet access and might be time-consuming. 

docker build -t hiclib-notebook .
