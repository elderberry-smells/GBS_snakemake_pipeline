# Pre-demultiplexed workflow for paired end GBS

Use this process if you are using ALREADY demultiplexed samples.  
It will run through the same process as the original snakefile (trim, align, sort, QC, cleanup) but with fewer inputs to complete. 
It will remove all temporary files, and zip up your demultiplexed folder at the end.  You should be left with a 
directory of sorted BAMs to use in your next process.

### What you will need to run this pipeline
1. Demutliplexed samples inside a directory (all together).  The format will look something like this (note that the naming of the files needs to be .1.fastq and .2.fastq)
```
project_directory/    
|
└─── demultiplexed_fastq/
    |    sample1.1.fastq
    |    sample1.2.fastq
    |    ...
    |    sample100.1.fastq
    |    sample100.2.fastq
```

2. A reference genome (indexed)

Thats it.  

### How to run the Snakefile-demux
1. Navigate to the project directory containing the demultiplexed directory.
```
$ cd /home/username/path/to/project_directory/
```
3. Open the snakemake config file to edit it (2 lines)
```
$ nano ~/gbs/GBS_Snakemake_pipeline/config/config-demux.yaml
```
3. With <ins>absolute paths</ins>, put in the full path to your demultiplexed folder and reference genome
```
demux_dir: "/home/username/path/to/prjoect_directory/demultiplexed_fastq"
reference: "/home/username/reference/reference1/reference1.fasta"
```
4. Run the Snakefile-demux (2 ways)

    a. On the cluster  (set cores by changing pre_demux.sh file parameters)
  
    ```
    $ qsub ~/gbs/GBS_Snakemake_pipeline/workflow/pre_demux.sh
    ```
    
    b. On your own machine (set cores after -j)
    
    ```
    $ conda activate gbs
    $ snakmake -s ~/gbs/GBS_Snakemake_pipeline/workflow/Snakefile-demux -j 8
    ```
    
5. Done

###  Outputs
Here is what your directory will look like once the run is completed
```
project_directory/    
|
└─── log/
|    | bwa_map.log
|    | novosort.log
|    | trimmomatic.output.log
|    └─── multiqc_report/
|        | multiqc_report.html
|        └─── multiqc_data/
|
└─── mapped_reads/
|
└─── sorted_bams/
|    | sample1.bam
|    | sample1.bam.bai
|    | ...
|    | sample100.bam
|    | sample100.bam.bai
|
└─── trimmomatic/
|
| config_demux.yaml 
| demultixed_fastq.tar.gz
```
