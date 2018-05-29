FROM centos:7

RUN yum update -y && yum install -y \
    environment-modules \
    git \
    glibc-static \
    gcc \
    gcc-c++ \
    boost \
    boost-devel \
    cmake \
    make \
    zlib-devel \
    ncurses-devel \
    java-1.8.0-openjdk \
    java-1.8.0-openjdk-devel \
    unzip \
    bzip2 \
    gcc-gfortran \
    perl-List-MoreUtils \
    perl-IO-Compress \
    libcurl-devel \
    libxml2-devel \
    openssl-devel \
    gd-devel \
    cairo-devel \
    which

WORKDIR /tmp

RUN curl -O https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm && \
    curl -O http://linuxsoft.cern.ch/cern/centos/7/cern/x86_64/Packages/parallel-20150522-1.el7.cern.noarch.rpm && \
    yum -y --nogpgcheck localinstall *.rpm

RUN git clone https://github.com/BenLangmead/bowtie.git && cd bowtie && git checkout tags/v1.1.2 && make && make install

RUN git clone https://github.com/BenLangmead/bowtie2.git && cd bowtie2 && git checkout tags/v2.2.6 && make && make install

RUN mkdir -p /software/tophat/2.1.0 && \
    curl -O http://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.0.Linux_x86_64.tar.gz && \
    tar -xzf tophat-2.1.0.Linux_x86_64.tar.gz -C /software/tophat/2.1.0 --strip-components=1

RUN git clone https://github.com/samtools/htslib.git && cd htslib && git checkout tags/1.2 && cd .. && \
    git clone https://github.com/samtools/samtools.git && cd samtools && git checkout tags/1.2 && make && make install

RUN git clone https://github.com/lh3/bwa.git && cd bwa && git checkout tags/0.7.12 && make && cp bwa /usr/local/bin

RUN mkdir -p /software/picard/1.131 && \
    curl -L -O https://github.com/broadinstitute/picard/releases/download/1.131/picard-tools-1.131.zip && \
    unzip picard-tools-1.131.zip && mv picard-tools-1.131 /software/picard/1.131

RUN mkdir -p /software/R/3.2.2 && \
    curl -O https://cran.cnr.berkeley.edu/src/base/R-3/R-3.2.2.tar.gz && \
    tar -xzf R-3.2.2.tar.gz && cd R-3.2.2 && ./configure --prefix=/software/R/3.2.2 --without-readline --without-x && \
    make && touch doc/NEWS.pdf && make install

RUN mkdir -p /software/openmpi/1.6.3 && \
    curl -O https://www.open-mpi.org/software/ompi/v1.6/downloads/openmpi-1.6.3.tar.gz && \
    tar -xzf openmpi-1.6.3.tar.gz && cd openmpi-1.6.3 && ./configure --prefix=/software/openmpi/1.6.3 && \
    make && make install

RUN git clone https://github.com/shenlab-sinai/ngsplot.git && cd ngsplot && git checkout tags/2.47 && cp -r bin/* /usr/local/bin

RUN git clone https://github.com/alexdobin/STAR.git && cd STAR && git checkout tags/2.5.2b && cp bin/Linux_x86_64_static/* /usr/local/bin

COPY docker/resources/modulefiles/ /etc/modulefiles/

COPY docker/resources/rsetup.R . 

RUN source /etc/profile.d/modules.sh && module load R && Rscript rsetup.R

COPY docker/resources/environment.yml .

RUN mkdir -p /software/anaconda2 && curl -L -O https://repo.continuum.io/archive/Anaconda2-2.4.1-Linux-x86_64.sh && \
    chmod 755 Anaconda2-2.4.1-Linux-x86_64.sh && bash Anaconda2-2.4.1-Linux-x86_64.sh -p /software/anaconda2 -f -b && \
    /software/anaconda2/bin/conda env update -f environment.yml --name root

# have to manually install bzip2 b/c the bioconda version of pysam depends on it
RUN source /etc/profile.d/modules.sh && module load python/anaconda && conda install -c conda-forge -y bzip2

RUN rm -rf *

RUN mkdir -p /projects/p20742/tools/bin

WORKDIR /projects/p20742/tools/bin

# copy all R and perl scripts into image
COPY bin/* ./

# install the stuff that the pipeline runs directly inside /projects/p20742
RUN curl -O http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed && chmod 755 bedToBigBed

RUN curl -O https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && unzip fastqc_v0.11.5.zip && chmod 755 FastQC/fastqc

RUN mkdir GATK_v3.6 && cd GATK_v3.6 && \
	curl -o GenomeAnalysisTK-3.6-0-g89b7209.tar.bz2 "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.6-0-g89b7209" && \
	bunzip2 GenomeAnalysisTK-3.6-0-g89b7209.tar.bz2 && tar xvf GenomeAnalysisTK-3.6-0-g89b7209.tar && \
	mv resources/* .

COPY docker/resources/SICER_V1.1.tgz .

RUN tar xvfz SICER_V1.1.tgz && cd SICER_V1.1/SICER && \
	find . -name '*.sh' -print | xargs sed -i 's|/home/data/SICER1.1|/projects/p20742/tools/SICER_V1.1|g'

RUN curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip && unzip Trimmomatic-0.33.zip

RUN git clone https://github.com/taoliu/MACS.git MACS-1.4.2 && cd MACS-1.4.2 && git checkout tags/v1.4.2 && python setup.py install --prefix /projects/p20742/tools/MACS-1.4.2 && ln -s /projects/p20742/tools/MACS-1.4.2/lib/python2.7 /projects/p20742/tools/MACS-1.4.2/lib/python2.6

RUN git clone https://github.com/deweylab/RSEM.git RSEM-1.3.0 && cd RSEM-1.3.0 && git checkout tags/v1.3.0 && make && \
    source /etc/profile.d/modules.sh && module load R && make ebseq

RUN curl -O https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_v0.11.4.tar.gz && tar -zxf fastq_screen_v0.11.4.tar.gz

COPY docker/resources/fastq_screen.allRefs.conf fastq_screen_v0.11.4/
# have to symlink perl executable in order for buildPipelineScripts.pl to run
RUN mkdir -p /software/activeperl/5.16/bin && ln -s /usr/bin/perl /software/activeperl/5.16/bin/perl

WORKDIR /projects/p20742/tools

COPY docker/resources/wrapper.sh .

RUN chmod 755 wrapper.sh

RUN ln -s bin/buildPipelineScripts.pl buildPipelineScripts.pl

ENTRYPOINT ["/projects/p20742/tools/wrapper.sh"]
