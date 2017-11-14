FROM centos:7

RUN yum update -y && yum install -y \
    environment-modules \
    git \
    glibc-static \
    gcc \
    gcc-c++ \
    boost \
    cmake \
    make

RUN curl -o /tmp/bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm && \
    yum -y --nogpgcheck localinstall /tmp/bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm && \
    rm /tmp/bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm

WORKDIR /tmp

RUN git clone https://github.com/BenLangmead/bowtie2.git && cd bowtie2 && git checkout tags/v2.2.6 && make && make install
