# Running the Ceto pipeline via Docker

The Ceto pipeline was originally developed to run on Northwestern University's
Quest HPC cluster. However, it can be run anywhere with via a Docker container.
This repository contains a `Dockerfile` that can be used to create a Docker
image and run the pipeline.

If you need to install Docker, documentation can be found at
[https://docs.docker.com/](https://docs.docker.com/).

## Building the image

Assuming you are working on a host with Docker already installed and running,
you can build the image with this command:

	docker build -t ceto:1.0 .


## Downloading the reference data

In order to align your sample data, you will need to download the reference data
that Ceto needs to a directory where it can be mounted by the Docker container.
The data is available in a public S3 bucket and in the example below we create a
new directory `/anno` and download it there:

	mkdir /anno
	aws s3 cp s3://ngsbartom-ceto-data/anno /anno --recursive

**NB:** If you are running this pipeline on an AWS EC2 instance, make sure that
the IAM Instance Profile of your instance has at least the
AmazonS3ReadOnlyAccess managed policy attached to it.

## Running the pipeline in a container

The pipeline has significant computing resource requirements. Ideally, it should
be run on a machine with at least as many cores as you have samples to align, so
that each alignment can run on its own core. The minimum recommended AWS EC2
instance size to run the pipeline on is a t2.2xlarge (8 cores, 32GB RAM).
Additionally, the reference data requires about 150GB of hard drive space and
another 50GB of scratch space is recommended.

1. First upload your sample data to a folder on the host where you will run your
   container, such as `/data`. In there, place the comparisons file (e.g.
   `comparisons.csv`)

2. Then create an output directory, e.g.

		mkdir /output

3. Finally, run the docker container. Arguments to the `docker run` command will
   be passed through to the `buildPipelineScripts.pl` script that runs the
   pipeline.
   
   **NB: Only the `build*` arguments to `buildPipelineScripts.pl` are supported,
   e.g. `buildAlign 1` and `buildPeakCaller 1`. The `run*` arguments are NOT
   supported, as they require the presence of an HPC scheduler. Instead, the
   scripts created by the pipeline are executed automatically in the proper
   order.**
      
   **NB: The volume mount path for the reference data inside the container
    must be `/projects/p20742/anno` in order for the pipeline to run.**

    **NB: The volume mount paths for the data and output directories must be
    the same both inside the container and on the host, as otherwise symlinks
    created by the pipeline will not work.**

	### RNA Example:

		docker run -it --rm -v /anno:/projects/p20742/anno -v /data:/data -v /output:/output ceto:1.0 \
		-t RNA -o /output -g mm10 -f /data -c /data/comparisons.csv -buildAlign 1 -buildEdgeR 1

	### ChIPseq example:

		docker run -it --rm -v /anno:/projects/p20742/anno -v /data:/data -v /output:/output ceto:1.0 \
		-t chipseq -o /output -g sacCer3 -f /data \
		-chip /data/<sample.csv> -buildAlign 1 -buildPeakCaller 1

When the run completes, the output can be found in the `/output` directory.

## Benchmarks

In January 2018, using the RNA command line above and the sample referenced in
[RNAseqAnalysisWithCETO.txt](RNAseqAnalysisWithCETO.txt), the pipeline completed
in the following times:

| Instance Type | Time (minutes) | Est. Cost* |
| ------------- | -------------- | ---------- |
| m4.2xlarge    | 465            | $3.10      |
| c4.4xlarge    | 256            | $3.40      |
| c4.8xlarge    | 158            | $4.19      |

The pipeline did not finish on a c4.2xlarge instance.

*EC2 on demand instance time only, us-east-2 region 1/4/2018.

## FAQs

**I get an Access Denied error when trying to download the reference data from
S3!**

This is probably because you either don't have your AWS CLI credentials set up
properly, or if you're running the pipeline on an EC2 instance, the instance's
IAM profile may not have the AmazonS3ReadOnlyAccess policy attached.

