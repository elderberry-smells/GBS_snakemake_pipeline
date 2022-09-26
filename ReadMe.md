<!-- [![AAFC](https://avatars1.githubusercontent.com/u/4284691?v=3&s=200)](http://www.agr.gc.ca/eng/home/)-->

# Paired End Genotype By Sequencing (GBS) Pipeline

> A bioinformatics package that includes worflows to generate bams and a vcf file from multiplexed GBS paired end fastq files.
  
> This tool utilizes the workflow management system [Snakemake](https://snakemake.readthedocs.io/en/stable/) and allows you to accomplish a complete pipeline process with just a few command line calls, even on a computing cluster.

An Example Workflow:

![workflow](workflow/resources/images/BAM_workflow.jpg?raw=true "Workflow")

## Outline

- Generate GBS data from Paired-End Fastq's generated by Illumina sequencing technologies
- Split multiplexed libraries into separate fastq's for each sample (read 1 and read 2) for in-line barcoded libraries.
- Process GBS libraries that are made of different species, or samples that require different reference genomes.
- End to end processing with only 2 command line calls to GridEngine. One for producing BAMs, one for producing VCF's
- Quality control - statictics produced for each sample and summarized in a MultiQC report 
- Designed to run on high performance clusters
- Clean data management.  This pipeline will remove all the intermediate files, and zip up the demultiplex files once complete to save disc space on your system.
- Custom iGenomics Riptide GBS library processing
- Ease of use - you need to only provide the tool with paths to a sample sheet, barcode file and a reference genome


## Table of Contents
- [Quick Use Guide](#quick-use-guide)
	- [Generating Bams](#generating-bams)
	- [Generating VCF](#generating-vcf)
- [Installation](#installation)
- [Supporting Files](#supporting-files)
- [Features](#features)
- [Team](#team)
- [License](#license)

## Quick Use Guide
<ins>ATTENTION</ins> Already have demultiplexed fastq's?  Try using the [Demux Pipeline](workflow/resources/ReadMes/pre-demux_workflow.md)

This quick use guide is assuming you have followed the installation.  If not, proceed first to [installations](#installation).

Once your environment is set up you can start running the pipeline!

Start at the beginning with Generating BAMs, or if you have BAMs you want to process, move on to generating VCFs!



#### Generating BAMs
Take paired end Fastq files (R1 and R2) and produce sorted BAM files and BAM index files
1. Create and save sample sheet into cluster - generally helpful to put it in same directory as fastq files.  For format [see samplesheet](workflow/resources/ReadMes/Supporting_Files.md) in supporting file readme.
2. Change directory to where fastq files are located
    
    ```$ cd directory/to/samples.fastq.gz```
3. Update the config/config.yaml file. In terminal, input the following commands. 
   ```
   $  nano ~/gbs/GBS_snakemake_pipeline/config/config.txt
   ```
   manually change the fields for your project
   ```
   sample_read1: "/input/absolute/path/to/sample_R1.fastq.gz"
   samplesheet: "/input/absolute/path/to/samplesheet.txt or samplesheet.csv"
   barcodefile:  "input/path/to/workflow/resources/barcodes/barcodes.txt"
   ```
    `cntrl-o` to save.  (press enter)

    `cntrl-x` to exit nano editor

4. Run the pipeline by submitting the process via shell script to the GridEngine queue in the terminal.

    ```shell script
    $ qsub ~/gbs/GBS_snakemake_pipeline/workflow/gbs.sh
    ```
5.  Done!   

*Notes*:

- This will give you a job number if submitted properly.  
- The job will produce snakemake stats (which rule is running) in your current working directory and be named `gbs.sh.e(jobnumber)`
- You can check to make sure the program is running by typing `qstat -f`.  If not, check the `gbs.sh.e(jobnumber)` for why.
- Runtime ~ 1 to 2 days for a standard output from Illumina HiSeq.  
- Output of the tool will be in the directory where your fastq files were located. 

#### Generating VCF
Take a list of BAM directories, and produce a single or multiple VCF files from all the BAM files included in those directories.
1. Update the config/config_vcf.yaml file for the 4 fields.
- For format [see config files](workflow/resources/ReadMes/Supporting_Files.md) in supporting file readme.
- *note* if you want to produce a seperate VCF for each chromosome in the species, put True into split: field.  This will speed up the process as each VCF will be written in parallel.  If your reference has scaffolds, it will produce a seperate VCF just for scaffolds (all together)

    ```shell script
    $  nano ~/gbs/GBS_snakemake_pipeline/config/config_vcf.txt
    ```
    ```
    sample_directories:
      - "/absolute/path/to/some/sorted_bams"
      - "/absolute/path/to/some/more/sorted_bams"
      - "/absolute/path/to/even/more/sorted_bams"
    
    # produce a single whole genome VCF (False) or split into 1 VCF for each chromosome (True) -- case sensitive
    split: True

    outfile: "/absolute/path/to/output/vcf/directory"
    
    reference: "/absolute/path/to/reference_genome.fasta"
    ```
    `cntrl-o` to save
    
    `cntrl-x` to exit
    
2. (optional) Create a directory for your VCF output (the same directory name as the config-vcf.yaml output path) and cd into that directory so your snakemake ouput log lives with the VCF file.  The name of the VCF will be whatever you name the output folder!
    
    ```shell script
    $ mkdir vcf_output/ && cd vcf_output
    ```
    
2. Run the pipeline by submitting the process via shell script to the GridEngine queue in the terminal.

    ```shell script
    $ qsub ~/gbs/GBS_snakemake_pipeline/workflow/vcf.sh
    ```

3. Done!

(optional) 4. Process the VCF file 
Check out the [VCF_Analysis_tools ReadMe](https://github.com/elderberry-smells/GBS_snakemake_pipeline/blob/master/workflow/resources/ReadMes/VCF_Analysis_Tools.md) to see a couple tools to process and visualize your VCF
- A filtering program detailed in the the following ReadMe can be used to filter the data for a variety of parameters
- Determine coverage of SNPs across the genome and graph that using the windowed_variance.py script


*Notes*:
- This will give you a job number if submitted properly.
- You can check to make sure the program is running by typing `qstat`.  Look for `vcf.sh` in the list.  If not, check the `vcf.sh.e(jobnumber)` for why.
- Final outputs will be in the folder you designated in the config_vcf.yaml file


## Installation

1. In your home directory on the cluster/local computer create a directory for the pipeline to live called `gbs` 
    ```shell script
    $ mkdir gbs
    $ cd gbs/
    ```

2. Clone the repository to that gbs folder you created on the cluster 
- You will need git installed on your cluster account to accomplish this.  Can be done using `conda install -c anaconda git` on your cluster account.

    ```shell script
    $ git clone https://github.com/elderberry-smells/GBS_snakemake_pipeline.git
    ```
    press `enter`, your local clone should be created with the following pathways for the pipeline:
    ```
    ~/gbs/GBS_snakemake_pipeline/snakefiles
    ~/gbs/GBS_snakemake_pipeline/config/
    ~/gbs/GBS_snakemake_pipeline/workflow/envs/    
    ~/gbs/GBS_snakemake_pipeline/workflow/rules/
    ~/gbs/GBS_snakemake_pipeline/workflow/scripts/
    ~/gbs/GBS_snakemake_pipeline/workflow/resources/
    ```
3. Install the GBS snakemake environment on your computer
- Create an environment for snakemake using the `workflow/envs/gbs-spec-file.txt` file from the repo you downloaded, then activate it.
    ```shell script
    $ conda create --name gbs --file ~/gbs/GBS_snakemake_pipeline/workflow/envs/gbs-spec-file.txt
    $ conda activate gbs
    ```
- When that is done running, you should see `(gbs)` on the left side of your command line
4a. Install Novocraft 3.10 if you don't have it (https://www.novocraft.com/support/download/) 
- Untar with command tar -xzf novo…tar.gz
4b.  Allow Novosort to be a callable program within gbs environment
- move that entire novocraft directory into the recently created gbs environment, and add that path to your bashrc profile 
    ```shell script
    $ mv ~/gbs/GBS_Snakemake_pipeline/workflow/resources/novocraft ~/miniconda3/envs/gbs/bin
    $ nano ~/.bashrc
    ```
- at the bottom of that file, add the following line:
```shell script
export PATH="~/miniconda3/envs/gbs/bin/novocraft":$PATH
```
- cntrl + O to save (enter)
- cntrl + X to close

- make the novosort file executable by changing permissions
    ```shell script
    $ chmod 777 ~/miniconda3/envs/gbs/bin/novocraft/novosort
    ```

- everything should be installed to run the pipeline now.
     
#### bioinformatics dependencies
all of these dependencies should be now be installed and callable by typing in their names in terminal

`python=3.6`
`snakemake=4.0`
`trimmomatic=0.39`
`bwa=0.7.17`
`samtools=1.9`
`bcftools=1.8`
`Novosort=2.00.00`

### Testing to see if the pipeline is working
> If you have data and a samplesheet to run, test to see if the snakemake process is working

Update the `config.yaml` file to direct the pipeline to the resources it needs.
```
$ nano ~/gbs/GBS_snakemake_pipeline/config/config.yaml
```

Invoke the snakemake process in test mode to see if outputs are being called correctly.
```shell script
$ snakemake -s ~/gbs/GBS_snakemake_pipeline/Snakefile -np
```
the snakemake process should show you all it intends to do with the samples located in the samplesheet, and how many process for each module/rule will run.  

example:

```
Job counts:
	count	jobs
	1	all
	184	bwa_map
	1	demultiplex
	184	novosort
	184	trimmomatic
        1       fastqc
        1       multiqc
	1	zipup_demux
	1	zipup_unmatched
	557
```
#### Supporting Files
A list of user inputs to run the program can be found in the [Supporting Files](workflow/resources/ReadMes/Supporting_Files.md) ReadMe.

This includes the config files, samplesheets, and barcode files.  

Refer to the docs for formatting of the files.  

Examples for samplesheets can be found in `workflow/resources/example_samplesheet.csv`

#### Features
Brief documentation on each step in this pipeline can be found in the [Features ReadMe](workflow/resources/ReadMes/Features.md)

The protocol for library construction is also included in the [resources folder](workflow/resources/oligos)

# Team
- Author:  Brian James (brian.james4@canada.ca)
- Testing and Improvements:  Dr. Jana Ebersbach (jana.ebersbach@canada.ca)
- Lead Scientist: Dr. Isobel Parkin (isobel.parkin@canada.ca)


## License
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
