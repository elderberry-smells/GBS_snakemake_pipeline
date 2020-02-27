<!-- [![AAFC](https://avatars1.githubusercontent.com/u/4284691?v=3&s=200)](http://www.agr.gc.ca/eng/home/)-->

# Paired End GBS Pipeline

> A bioinformatics tool to generate a vcf file from multiplexed GBS paired end fastq files

> Utilizes snakemake to complete the process with one command line call

**maybe put in a graphic here (gif) to record how the tools works?  Or put in the dag plot from snakemake**

![workflow](https://github.com/elderberry-smells/GBS_snakemake_pipeline/tree/master/workflow/resources/images/pipeline_workflow_DAG.jpg)

## Table of Contents

- [Installation](#installation)
- [Quick Use Guide](#quick-use)
- [Features](#features)
- [Usage](#usage)
- [running the pipeline](#Executing)
- [Team](#team)
- [License](#license)

## Installation
> download the repository using [git clone](#clone) or by clicking `clone or download` on the main page of the repository
and dowloading the zip file.  
- in your home directory on the cluster/local computer create a directory for the pipeline to live called `gbs` 
```shell script
$ mkdir ~/gbs & cd ~/gbs/
```
### clone 
> you will need git installed on your computer to accomplish this.  Can be done using `conda install -c anaconda git` on your cluster account.
- Use the `git clone` process to make an exact copy of the current repository
- clone this repository to your local machine using the following:

```shell script
$ git clone https://github.com/elderberry-smells/GBS_snakemake_pipeline.git
```
- press `enter`, your local clone should be created with the following pathways for the pipeline:

    - `~/gbs/GBS_snakemake_pipeline/snakefile`
    - `~/gbs/GBS_snakemake_pipeline/config/`
    - `~/gbs/GBS_snakemake_pipeline/workflow/envs/`    
    - `~/gbs/GBS_snakemake_pipeline/workflow/rules/`
    - `~/gbs/GBS_snakemake_pipeline/workflow/scripts/`
 
### installing the GBS snakemake environment on your computer
- ensure you are running the most recent version of conda (if you can)
```shell script
$ conda update conda
``` 
- Create an environment for snakemake using the `workflow/envs/gbs.yaml` file from the repo you downloaded
```shell script
$ conda env create -f ~/gbs/GBS_snakemake_pipeline/workflow/envs/gbs.yaml
```
#### bioinformatics dependencies
- all of these dependencies should be installed from the environment being created except for novosort

`python=3.6`
`trimmomatic=0.39`
`bwa=0.7.17`
`samtools=1.9`
`bcftools=1.8`
`Novosort=2.00.00`

- move the novosort folder to your newly created gbs environment bin
    
    ```shell script
    $ mv GBS_snakemake_pipeline/workflow/resources/novocraft/ ~/miniconda3/envs/gbs/bin/
    ```
      
    - add that folder to your profile so novosort is a callable command
    
    ```shell script
    $ nano ~/.bashrc 
    ```
    - add the path to the novocraft folder to the bottom of the file, save, and exit.  Change the path below to the correct one on your machine
    
    `export PATH="$PATH:~/miniconda3/envs/gbs/bin/novocraft"` 

### Some minor changes to the snakefile and rule files to include your user name
- You will have to edit the file to reflect your username in the scripts before it will run.  You can do this with whichever editor you would like, `nano` would be fine for these small edits.  
- replace /home/AAFC-AAC/`your_user_name`/gbs... with your actual username
##### Snakefile
-  line 30
-  line 52 - 56 
##### workflow/rules/demultiplex.smk
- line 10
##### workflow/rules/trimmomatic.smk
- line 5

## Quick Use Guide

> Once your environment is set up, and the pathways reflect your user name, you can start running the pipeline.

- input required:
	- `sample_R1.fastq.gz` `sample_R2.fastq.gz` (no need to unzip)
	- `samplesheet.txt` for format see samplesheet below

- update the config/config.yaml file.  Use absolute paths.
```shell script
$  cd gbs/GBS_snakemake_pipeline
$  nano config/config.txt
```
`cntrl-o to save [enter]`
`cntrl-x to exit`

- run the pipeline by submitting the process to the gridengine queue

```shell script
$ qsub workflow/resources/gbs.sh
```

## Features

This tool utilizes Snakemake to create a cradle-to-grave GBS analysis for paired end reads from Illumina sequencing platforms (2x150).

The tool will commit the following steps in the pipeline, all of which are modular additions to the snakefile and can be swapped if needed with other tools.

#### demultiplex
In paired end reads, the Fastq read 1 houses the unique identifier barcodes for demultiplexing.  The barcodes in this tool are designed based on the 
[Poland et. al 2012 GBS protocol](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0032253), and the barcodes being used (up to 384 unique barcodes for multiplexing) can be found in `workflow/resources/barcodes_384.txt`

the barcodes themselves are anywhere between 5 and 14 bp in length, and end in `TCGA` signalling the start to the 
actual sample sequence.  This tool will demultiplex read 1 and read 2 into seperate files `sample_id.1.fq` and `sample_id.2.fq`

#### trimmomatic
The next step is to pipe the demultiplexed files into the trimmomatic tool.  
> Trimmomatic is a fast, multithreaded command line tool that can be used to trim and crop
Illumina (FASTQ) data as well as to remove adapters. These adapters can pose a real problem
depending on the library preparation and downstream application.

outputs of this file will be `sample_id.1.paired` `sample_id.1.unpaired` `sample_id.2.paired` `sample_id.2.unpaired`

#### Alignment and Sorting

trim files are passed through bwa mem for alignment to the reference genome, and then sorted and indexed using novosort

#### Generating SNP calls and VCF

calls are generated using samtools mpileup, and visualized in a VCF file using bcftools 

## Usage

The snakemake tool is structured in a way that you need only activate the environment and run the snakefile  after you update the `config/config.yaml` to direct the program.

The config file is 4 lines, all required inputs for the snakefile to work properly

```
sample_read1: "/absolute/path/to/sample_R1.fastq"
samplesheet: "/absolute/path/to/samplesheet.txt"
barcodefile:  "workflow/resources/barcodes/barcodes.txt"
reference_file: "/absolute/path/to/reference.fasta
```

The sample sheet is a tab delimited txt file with 3 columns.  An example of a samplesheet is shown below:

```
Sample_number	Index_name	Sample_ID
1	gbsx001	Sample1	
2	gbsx002	Sample2	
3	gbsx003	Sample3	
4	gbsx004	Sample4	
5	gbsx005	Sample5	
6	gbsx006	Sample6	
7	gbsx007	Sample7	
8	gbsx008	Sample8	
9	gbsx009	Sample9
10	gbsx010	Sample10
```


The barcode file, as seen in the barcode folder `workflow/resources/barcodes/` are 2 columns, tab delimited txt files.  This file does not have any headers

`col1: index_name` `col2: barcode`

example barcode file:

```
gbsx001	TGACGCCATGCA
gbsx002	CAGATATGCA
gbsx003	GAAGTGTGCA
gbsx004	TAGCGGATTGCA
gbsx005	TATTCGCATTGCA
gbsx006	ATAGATTGCA
gbsx007	CCGAACATGCA
gbsx008	GGAAGACATTGCA
gbsx009	GGCTTATGCA
gbsx010	AACGCACATTTGCA
```

## Executing

### Using interactive node (not recommended for large files)

Update the `config/config.yaml` file to direct the pipeline to the resources it needs.

If running on interactive node, qlogin with requested cores before activating environment 

`qlogin -pe smp 16`

then activate the GBS pipeline and change directory to where the snakefile is located  
```shell script
$ conda activate gbs
$ cd ~/gbs/GBS_snakemake_pipeline/
```

You can invoke the process by running one command, and if you have the ability to request cores using the -j option

- test to see if the snakefile is reading your data correctly

```shell script
$ snakemake -np
```
the snakemake process should show you all it intends to do with the samples locate din the samplesheet, and how many process for each module.  example:

```
Job counts:
	count	jobs
	1	all
	184	bwa_map
	1	demultiplex
	184	novosort
	1	samtools_call
	184	trimmomatic
	555
```

- run the program with requested cores on local machine or interactive node on cluster (use nohup in case you need to close terminal) 
```shell script
$ nohup snakemake -j 16
``` 

###  Running the program on the cluster with qsub using shell script from `workflow/resources/gbs.sh`
- from the head node, navigate to the folder with the snakefile.  Make sure your `config/config.yaml` is updated first!
```shell script
$ cd gbs/GBS_snakemake_pipeline/
$ qsub workflow/resources/gbs.sh
```
- this will give you a job number if submitted properly
- you can check to make sure the program is running by typing `qstat -f`
- runtime = 5ish days.  Working on a faster demux with threading to cut a day off.  
  
