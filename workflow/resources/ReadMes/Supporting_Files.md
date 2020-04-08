<!-- [![AAFC](https://avatars1.githubusercontent.com/u/4284691?v=3&s=200)](http://www.agr.gc.ca/eng/home/)-->

# Pipeline Supporting Files 
There are a few different files used to support the program that are inputs by the user.  This documentation will detail the formats of each.

1. [Config Files](#Config-file)
2. [Sample sheet](#Samplesheet) 
2. [Barcode File](#Barcode-file)


## Config file
- all config files are found in the `config` directory in the pipeline
- Make sure you use absolute paths. 
- Your paths need to be covered by "" at either end (this denotes a string in python)
- You can't short-form the paths with ~ as this is a value in the program, not a variable

#### 1.  BAM Generation (config.yaml)

The config file is 3 lines, all required inputs for the snakefile to work properly. 
- We only need R1 fastq file in config.
- that means your paired end (R2 file) needs same name (with R2 replacing R1) for this to work properly.
- fastq files should remain gzipped!
```
sample_read1: "/absolute/path/to/sample_R1.fastq.gz"
samplesheet: "/absolute/path/to/samplesheet.txt"
barcodefile:  "absolute/path/to/barcodes_192.txt"
```

#### 2.  VCF Generation (config_vcf.yaml)

The config_vcf file has 3 fields to fill out, all are required for program.
- sample directories can be multiple lines, each pointing to a directory of bams (with correcsponding index files .bai)
- start each new line of sample directories with 2 spaces, a dash, a space, then path.
- outfile is the directory where you want the VCF to be generated.  If it doesn't exist, it will be created.
- reference should be directly to an unzipped reference.fasta file.  
- reference needs to be indexed using `bwa index reference.fast` if not already done!

```
sample_directories:
  - "/absolute/path/to/some/sorted_bams"
  - "/absolute/path/to/some/more/sorted_bams"
  - "/absolute/path/to/even/more/sorted_bams"
outfile: "/absolute/path/to/output/vcf/directory"
reference: "/absolute/path/to/reference_genome.fasta"
```

## Samplesheet
The sample sheet is a tab delimited txt file or a csv file with 4 columns.  

**Your sample sheet must be formatted in same way, so make sure you have a header and be wary of column order!**

The samples can have any number of different references, each sample can even have a different genome if required.  An example of a samplesheet with 2 references is shown below:


```
Sample_number	Index_name	Sample_ID   Reference_path
1	gbsx001	Sample1	    /absolute/path/to/reference1/reference_genome1.fasta
2	gbsx002	Sample2     /absolute/path/to/reference1/reference_genome1.fasta
3	gbsx003	Sample3	    /absolute/path/to/reference1/reference_genome1.fasta
4	gbsx004	Sample4	    /absolute/path/to/reference1/reference_genome1.fasta
5	gbsx005	Sample5	    /absolute/path/to/reference1/reference_genome1.fasta
6	gbsx006	Sample6	    /absolute/path/to/reference2/reference_genome2.fasta
7	gbsx007	Sample7	    /absolute/path/to/reference2/reference_genome2.fasta
8	gbsx008	Sample8	    /absolute/path/to/reference2/reference_genome2.fasta
9	gbsx009	Sample9     /absolute/path/to/reference2/reference_genome2.fasta
10	gbsx010	Sample10    /absolute/path/to/reference2/reference_genome2.fasta
```

example samplesheets (a csv and a txt version) can be found in the resources directory. 

## Barcode file
The barcode file is a tab delimited txt file, matching the index name to its corresponding barcode sequence.

Some barcode files are included in the download and can be found barcode folder `workflow/resources/barcodes/`

This file does not have any headers

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
