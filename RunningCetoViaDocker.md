# Running the Ceto pipeline via Docker

# Building the image
	`docker build -t nuitrcs/ceto:0.1 .`

# Downloading the reference data
	`aws s3 cp s3://ngsbartom-ceto-data/anno . --recursive`

# Running the pipeline in a container
	- create output folder, data folder.
		- put comparisons.csv in output folder

	- mounting the data
		rna
		`docker run -it --rm -v /mnt/anno:/projects/p20742/anno -v /mnt/fastq:/data -v /mnt/output:/output nuitrcs/ceto:0.1`

		chipseq
		`docker run -it --rm -v /mnt/anno:/projects/p20742/anno -v /mnt/testChIP:/data -v /mnt/testChIP-output:/output nuitrcs/ceto:0.1`

	- ONE STEP RUN:
		docker run -it --rm -v /mnt/anno:/projects/p20742/anno -v /mnt/testChIP:/data -v /mnt/testChIP-output:/output nuitrcs/ceto:0.1 -t chipseq -o /output -g sacCer3 -f /data/GSE53241/ -chip /data/chipseqDesc.GSE53241.csv -buildAlign 1 -buildPeakCaller 1

		When this completes, output will be in the /mnt/testChIP-output directory.

