FROM centos:7

RUN yum update -y && yum install -y \
    environment-modules \
    git \
    glibc-static \
    gcc \
    gcc-c++ \
    boost \
    cmake \
    make \
    zlib-devel \
    ncurses-devel \
    java-1.8.0-openjdk \
    java-1.8.0-openjdk-devel \
    unzip \
    bzip2 \
    gcc-gfortran

WORKDIR /tmp

RUN curl -O https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm && \
    yum -y --nogpgcheck localinstall bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm

RUN git clone https://github.com/BenLangmead/bowtie.git && cd bowtie && git checkout tags/v1.1.2 && make && make install

RUN git clone https://github.com/BenLangmead/bowtie2.git && cd bowtie2 && git checkout tags/v2.2.6 && make && make install

RUN mkdir -p /software/tophat/2.1.0 && \
    curl -O http://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.0.Linux_x86_64.tar.gz && \
    tar -xzf tophat-2.1.0.Linux_x86_64.tar.gz -C /software/tophat/2.1.0 --strip-components=1

RUN git clone https://github.com/samtools/htslib.git && cd htslib && git checkout tags/1.2 && cd .. && \
    git clone https://github.com/samtools/samtools.git && cd samtools && git checkout tags/1.2 && make && make install

RUN git clone https://github.com/arq5x/bedtools2.git && cd bedtools2 && git checkout tags/v2.18.0 && make && \
    mkdir -p /software/bedtools/2.18.0 && cp bin/* /software/bedtools/2.18.0

RUN git clone https://github.com/lh3/bwa.git && cd bwa && git checkout tags/0.7.12 && make && cp bwa /usr/local/bin

RUN mkdir -p /software/picard/1.131 && \
    curl -L -O https://github.com/broadinstitute/picard/releases/download/1.131/picard-tools-1.131.zip && \
    unzip picard-tools-1.131.zip && mv picard-tools-1.131 /software/picard/1.131

COPY resources/environment.yml .

RUN mkdir -p /software/anaconda2 && curl -L -O https://repo.continuum.io/archive/Anaconda2-2.4.1-Linux-x86_64.sh && \
    chmod 755 Anaconda2-2.4.1-Linux-x86_64.sh && bash Anaconda2-2.4.1-Linux-x86_64.sh -p /software/anaconda2 -f -b && \
    /software/anaconda2/bin/conda env update -f environment.yml --name root

RUN mkdir -p /software/R/3.2.2 && \
    curl -O https://cran.cnr.berkeley.edu/src/base/R-3/R-3.2.2.tar.gz && \
    tar -xzf R-3.2.2.tar.gz && cd R-3.2.2 && ./configure --prefix=/software/R/3.2.2 --without-readline --without-x && \
    make && touch doc/NEWS.pdf && make install

RUN rm -rf *
