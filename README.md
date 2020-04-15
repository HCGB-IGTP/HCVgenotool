# HCVgenotool instructions

## Software requeriments
	- Python 3 (tested on Python 3.6.3). Required python packages:
		* cutadapt (tested on 1.16)
		* numpy (tested on 1.14.1)
		* biopython (tested on 1.70)
		* pandas (tested on 0.22.0)
	- BWA (Burrows-Wheeler Aligner). Tested on 0.7.17-r1188. 
	  [Download here](https://sourceforge.net/projects/bio-bwa/files/)
	- ea-utils/fastq-join. Tested on 1.04.807.
	  [Download here](https://expressionanalysis.github.io/ea-utils/)
	- fastq\_quality\_filter from fastx-toolkit. Tested on 0.0.13 .
	  [Download here](http://hannonlab.cshl.edu/fastx_toolkit/download.html)

This tool was was designed to run from a linux terminal. Tested on Ubuntu 17.10 (x86-64 architecture).

## Run instructions
```
$cd HCV_genotyping/source
$./genotyping_pipeline.sh \
    ~/FASTQ \
    ~/primers/primers.fasta \
    ~/references/HCV_genotypes.fasta \
    ~/results
```

## Required arguments
	FASTQDIR: Folder where all the fastq files are contained. Only Illumina
 		pair-end sequencing from PCR amplicons is accepted. One sample pair per
 		file.
 		IMPORTANT: *.fastq or *.fastq.gz are accepted. Read 1 must be specified
		as *R1* and read 2 as *R2* in the filename. Sample name will be taken
 		from the string before the first "\_" character that should be present in
 		the filename.
 		Examples of valid fastq filenames:
 			2\_S1\_L001\_R1\_001.fastq.gz -> sample name: "2", read 1.
 			sample1\_R2.fastq          -> sample name: "sample1", read 2.
 	ADAPTERS\_FILE: path to a fasta file containing adapter sequences used to
 		create amplicons for sequencing. Only two adapters allowed, named as
 		"\>Forward" and "\>Reverse" ("\>F" and "\>R" or "\>f" and "r" is sufficient).
 		Adapter sequence has to span only one line.
 	REFERENCE: path to a fasta file containing HCV reference sequences.
		Sequence names should follow the convention TS\_XXXXX, where T is HCV
 		type, S is HCV subtype and XXXXX is the sequence identifier.
		HCV reference genomes are included in HCV\_references/HCV\_genotypes.fasta
		with date January 1st 2017, downloaded from [here](https://talk.ictvonline.org/ictv_wikis/flaviviridae/hepacivirus/m/hepacivirus-files/6789)
 	WORKINGDIR: path to a folder where output results. If WORKINGDIR does not
 		exist it will be created. If WORKINGDIR exists, program exits with
		ERROR message.

Results will be written in `~/results` folder (summary table at `~/results/results.csv`).

**IMPORTANT!:** please, always mantain the order in which arguments are given to `genotyping\_pipeline.sh` script.
