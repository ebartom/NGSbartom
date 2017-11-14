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
    java-1.8.0-openjdk

RUN curl -o /tmp/bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm && \
    yum -y --nogpgcheck localinstall /tmp/bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm && \
    rm /tmp/bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm

WORKDIR /tmp

RUN git clone https://github.com/BenLangmead/bowtie2.git && cd bowtie2 && git checkout tags/v2.2.6 && make && make install

RUN mkdir -p /software/tophat/2.1.0 && \
    curl -O http://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.0.Linux_x86_64.tar.gz && \
    tar -xzf tophat-2.1.0.Linux_x86_64.tar.gz -C /software/tophat/2.1.0 --strip-components=1

RUN git clone https://github.com/samtools/htslib.git && cd htslib && git checkout tags/1.2 && cd .. \
    && git clone https://github.com/samtools/samtools.git && cd samtools && git checkout tags/1.2 && make && make install
