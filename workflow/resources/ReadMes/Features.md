<!-- [![AAFC](https://avatars1.githubusercontent.com/u/4284691?v=3&s=200)](http://www.agr.gc.ca/eng/home/)-->

# Paired End Genotype By Sequencing (GBS) Pipeline - Features
- [split rule](#split)
- [demultiplex rule](#demultiplex)
- [Trim rule](#trimmomatic)
- [align rule](#alignment)
- [sort rule](#sorting)
- [Quality Control](#quality-control)
- [vcf generation](#generating-snp-calls-and-vcf-file)

### split
input:  

`fastq_R1.fastq.gz, fastq_R2.fastq.gz`

output: 
```
chunks/fastq_R1_aa.fq ... chunks/fastq_R1_zz.fq, 
chunks/fastq_R2_aa.fq ... chunks/fastq_R2_zz.fq
```

the program will count the lines in the sample_R1.fastq.gz file and divy that up into split files using `split` command (default 100M lines for split).  Currently can only accomplish split permuations from aa to zz.  Each file will come out of this stage into a newly created chunks/ folder.

### demultiplex
input: 
```
chunks/fastq_R1_aa.fq...chunks/fastq_R1_zz.fq, 
chunks/fastq_R2_aa.fq...chunks/fastq_R2_zz.fq
```

output: 
```
ref_name/demultiplex/sample1.1.fastq ... ref_name/demultiplex/sample384.1.fastq, 
ref_name/demultiplex/sample1.2.fastq ... ref_name/demultiplex/sample384.2.fastq
```

In paired end reads, the Fastq read 1 houses the unique identifier barcodes for demultiplexing.  The barcodes in this tool are designed based on the 
[Poland et. al 2012 GBS protocol](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0032253), and the barcodes being used (up to 384 unique barcodes for multiplexing) can be found in `workflow/resources/barcodes_384.txt`

This tool will demultiplex read 1 and read 2 into seperate files by matching the barcodes in the sequence and appending to the newly created sample files in chunks of 25 million sequences at a time (done through the usage of multiprocessing).  The [demux script](https://github.com/elderberry-smells/GBS_snakemake_pipeline/blob/master/workflow/scripts/PE_fastq_demultiplex.py) reads both read 1 and read 2 concurrently, so matching of header information in fastq is critical for the process.

parameters for demultiplex script:

`python3 PE_fastq_demultiplex.py -f fastq_R1.fastq.gz -b barcode_file.txt -s samplesheet.txt`

### trimmomatic
input: 
```
demultiplex/sample1.1.fastq ... demultiplex/sample384.1.fastq,
demultiplex/sample1.2.fastq ... demultiplex/sample384.2.fastq
```

output:

```
trimmomatic/sample_id.1.paired, trimmomatic/sample_id.1.unpaired
trimmomatic/sample_id.2.paired, trimmomatic/sample_id.2.unpaired
```

The next step is to pipe the demultiplexed files into the trimmomatic tool.  Trimmomatic is freely available and is a fast, multithreaded command line tool that can be used to trim and crop Illumina (FASTQ) data as well as to remove adapters. 

parameters used in trimmomatic:
```
trimmomatic PE -threads 16 -phred33 sample_id.1.fastq sample_id.2.fastq out.1.paired out.1.unpaired out.2.paired output.2.unpaired
ILLUMINACLIP:{trim_file}:2:40:15 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:55
```

### Alignment
input:
```
trimmomatic/sample_id.1.paired, trimmomatic/sample_id.1.unpaired,
trimmomatic/sample_id.2.paired, trimmomatic/sample_id.2.unpaired
```

output:

`mapped_reads/sample1.bam ... mapped_reads/sample384.bam`

The trimmed files are passed through the BWA MEM for alignment to the reference genome.  The reference genome is grabbed by looking into the original config file updated by the user.  Threading for this is set at 16 but can be upped if you reserved more cores by altering the workflow/resources/gbs.sh, as well as the threads: in the bwa_map.smk rule.

parameters used for bwa mem:

`bwa mem -t 16 reference.fasta trim1.fq trim2.fq | samtools view -Shbu > sample.bam`

### Sorting
input:

`mapped_reads/sample1.bam ... mapped_reads/sample384.bam`

output:
```
sorted_bams/sample1.bam ... sorted_bams/sample384.bam,
sorted_bams/sample1.bam.bai ... sorted_bams/sample384.bam.bai
```

This rule uses the freely available software Novosort (from Novocraft).  This will sort and index the bam filesfrom the bwa mem rule.  The sorted bam file will move into the SNP calling in the final rule

parameters used for novosort:

`novosort sample_id.bam --threads 16 --index --output sample_id.sorted.bam`

### Quality Control

#### fastqc

input:

`ref_dir/demultiplex/*.fastq`

output:
`ref_dir/log/fastqc/fastqc_report.html`

Run fastqc stats on all demultiplexed fastq files for each directory created in the pipeline

parameters:

`fastqc -o ref_dir/log/fastqc/ -t 16 ref_dir/demultiplex/*.fastq`

#### MultiQC

input:

`ref_dir/log/`

output:

`ref_dir/log/multiqc_report.html`

Scan through the log files from entire pipeline for each directory created and produce a summarized report for all samples.

parameters:

`multiqc ref_dir/log -o ref_dir/log/`

## Generating SNP calls and VCF file
The second part of the tool is the Snakefile-vcf script that creates a VCF file from a list of bam directories.

input:
- list of directories: `bam.list`
`input_directory_1/*.bam...input_directory_N/*.bam` 

output:

`output_directory/calls/prefix.vcf.gz`

Calls are generated using bcftools mpileup, and visualized in a VCF file using bcftools calls.  

parameters used in samtools mpileup:

`bcftools mpileup -a AD,DP -f reference.fasta -b bam.list | bcftools call -mv -Oz > output_directory/calls/prefix.vcf.gz`


