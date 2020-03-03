from string import ascii_lowercase
import os
import subprocess

def parse_samplesheet(sample_sheet):
    '''open samplesheet and grab relevant sample names.  Format of the sample sheet (tab delimited):
    Sample_number	Index_name	Sample_ID'''

    index = []
    samples = []
    with open(samplesheet, 'r') as samps:
        next(samps)  # skip the header line
        lines = samps.readlines()
        for line in lines:
            data = line.split()
            index.append(data[1])
            samples.append(data[2])
            
    return index, samples


def bam_list(path, sample_list):
    "make a list of bam files to run through the samtools mpileup"
    bamlist_name = f'{path}/bams.list'
    with open(bamlist_name, "w") as bam:
        for i in sample_list:
            line = f'{path}/sorted_reads/{i}.sorted.bam\n'
            bam.write(line)

    return bamlist_name


def generate_split_list(fastq_r1, fname, split_lines):
    "make a list of split files and their prefix for the split.smk rule"
    
    # generate a dictionary of all possible names in a prefix ex) aa, ab, ac...zx, zy, zz
    lets = {str(index): letter for index, letter in enumerate(ascii_lowercase, start=0)}
    
    total_names = []
    
    for i in range(26):  # first letter of prefix
        for x in range(26):  # second letter of prefix
            name = f'{fname}_{lets[str(i)]}{lets[str(x)]}.fq'
            total_names.append(name)
            
    # read the fastq_r1 file and determine how many lines there are
    wc_command = f'zcat {fastq_r1} | wc -l'
    line_count = float(subprocess.getoutput(wc_command))
    
    split_files = line_count/split_lines # total number of files that will be created
    
    name_list = total_names[:int(split_files)]  # slice only the names for the number of split files to be generated
    name_list2 = [split_fname.replace('_R1', '_R2') for split_fname in name_list]
    
    return name_list, name_list2
                      
    
# read in the config file to get pathways/samples into variables
configfile: "/home/AAFC-AAC/your_user_name/gbs/GBS_snakemake_pipeline/config/config.yaml"
read1 = config["sample_read1"]

# if the read number is in lowercase, make it uppercase
if '_r1' in read1:
    read1 = read1.replace('_r1', '_R1')
read2 = read1.replace('_R1', '_R2')  # make the read2 name (only difference is R1 vs R2)

samplesheet = config["samplesheet"]
barcode_file = config["barcodefile"]
ref_genome = config["reference_file"]

# get the pathway from the read1 fastq file, this is where we want to run the pipeline (set working dir)
path_name, base_name = os.path.split(read1)
file_handle = os.path.splitext(base_name)[0]
file_name = os.path.splitext(file_handle)[0]

# read the sample sheet and parse the data you need out of it into variables
index, samples= parse_samplesheet(samplesheet)

# how big do you want the chunk files to be.  default is 100M.  If memory is an issue, try smaller (like 10M).
split_number = 100000000  # must be a multiple of 4!!!!!

bams = bam_list(path_name, samples)  # create a list of bams to be run through mpileup at end stage

workdir: path_name

#######################################  RULES ##################################################

rule all:
    input:
        "calls/snps.raw.vcf.gz"  # final output is a single vcf file

include: "/home/AAFC-AAC/your_user_name/gbs/GBS_snakemake_pipeline/workflow/rules/split.smk"
include: "/home/AAFC-AAC/your_user_name/gbs/GBS_snakemake_pipeline/workflow/rules/demultiplex.smk"
include: "/home/AAFC-AAC/your_user_name/gbs/GBS_snakemake_pipeline/workflow/rules/trimmomatic.smk"
include: "/home/AAFC-AAC/your_user_name/gbs/GBS_snakemake_pipeline/workflow/rules/bwa_map.smk"
include: "/home/AAFC-AAC/your_user_name/gbs/GBS_snakemake_pipeline/workflow/rules/novosort.smk"
include: "/home/AAFC-AAC/your_user_name/gbs/GBS_snakemake_pipeline/workflow/rules/samtools_call.smk"
