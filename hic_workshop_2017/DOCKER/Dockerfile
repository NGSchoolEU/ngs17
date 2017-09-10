# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
FROM jupyter/minimal-notebook

USER root

# libav-tools for matplotlib anim
RUN apt-get update && \
    apt-get install -y --no-install-recommends libav-tools && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER $NB_USER

# Install Python 3 packages
# Remove pyqt and qt pulled in for matplotlib since we're only ever going to
# use notebook-friendly backends in these images
RUN conda config --add channels bioconda
RUN conda upgrade conda

RUN conda install --quiet --yes \
		'pigz=2.3' \
    'bowtie2=2.*' \
    'samtools=1.5*' \
		'pysam' \
    'pandas=0.19*' \
    'numexpr=2.6*' \
    'matplotlib=2.0*' \
    'scipy=0.19*' \
    'seaborn=0.7*' \
		'ipywidgets=6.0*' \
    'scikit-learn=0.18*' \
    'scikit-image=0.12*' \
    'cython=0.25*' \
    'statsmodels=0.8*' \
    'hdf5=1.8.17' \
    'h5py=2.6*' && \
    conda remove --quiet --yes --force qt pyqt && \
    conda clean -tips

# Adding some specific packages for Hi-C data analysis
RUN pip install https://bitbucket.org/mirnylab/mirnylib/get/tip.tar.gz
RUN	pip install https://bitbucket.org/mirnylab/hiclib/get/tip.tar.gz && \
	# Patch for fixing "long file name error"
	sed -i 's/".fastq.gz"/".gz"/g' /opt/conda/lib/python3.6/site-packages/hiclib/mapping.py

RUN cd /tmp && \
  wget --no-proxy https://github.com/nezar-compbio/lavaburst/archive/master.zip && \
  unzip master.zip && \
	cd lavaburst-master/ && \
  python setup.py install  --prefix='/opt/conda/'

# Adding ipython notebooks with code
#COPY ./*py* work/

# Activate ipywidgets extension in the environment that runs the notebook server
RUN jupyter nbextension enable --py widgetsnbextension --sys-prefix

# Import matplotlib the first time to build the font cache.
ENV XDG_CACHE_HOME /home/$NB_USER/.cache/
RUN MPLBACKEND=Agg python -c "import matplotlib.pyplot"

USER $NB_USER
