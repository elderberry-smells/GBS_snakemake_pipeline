from string import ascii_lowercase
import os
import subprocess
from pathlib import Path
import csv

#################################     Defined Functions      #########################################
def parse_samplesheet(sample_sheet):
    # determine the format of the samplesheet
    sample_path, sample_basename = os.path.split(sample_sheet)
    samp_ext = os.path.splitext(sample_basename)[1]

    try:
        samp_ext in ['.csv', '.txt', '.tsv']
    except NameError:
        print("Ooops!  Your file format for your samplesheet is not accepted for this pipeline.\n"
              "Change to one of the following formats:  .csv, .txt, .tsv")

    indexes = []
    samples = []
    ref_dirs = []
    refs = {}

    with open(sample_sheet, newline='') as samp_file:
        next(samp_file)
        if samp_ext == '.csv':
            reader = csv.reader(samp_file, delimiter=',')
        else:
            reader = csv.reader(samp_file, delimiter='\t')
        for row in reader:
            indexes.append(row[1])
            samples.append(row[2])
            ref_dirs.append(row[3].strip())
            # add to the reference dictionary:  format {sample_name: reference path}
            refs[row[2]] = row[3].strip()  # make sure there is no white space after line ends...

    return indexes, samples, ref_dirs, refs


def reference_samples(sample_list, reference_list):
    return {k: v for k, v in zip(sample_list, reference_list) if v != ''}


def output_lists(zipped_samples):
    """
    make a list of the outputs expected in the pipeline.
    :param zipped_samples: the samples and reference lists from config file zipped together
    :return: list for demux, trim, and bams, and sorted bams
    """
    demuxer = []
    bammer = []
    unique_ref_dirs = []
    for r, s in zipped_samples:
        if r == '':
            entry = f"unreferenced/demultiplex/{s}"
        else:
            direc = r.split('/')[-2]
            entry = f"{direc}/demultiplex/{s}"
            bam_entry = f"{direc}/mapped_reads/{s}"
            bammer.append(bam_entry)

            if direc not in unique_ref_dirs:
                unique_ref_dirs.append(direc)

        demuxer.append(entry)

    trimmer = [i.replace('demultiplex', 'trimmomatic') for i in demuxer]
    sorter = [i.replace('mapped_reads', 'sorted_bams') for i in bammer]


    return demuxer, trimmer, bammer, sorter, unique_ref_dirs


def generate_split_list(fastq_r1, fname, split_lines):
    """
    generate a dictionary of all possible names in a prefix ex) aa, ab, ac...zx, zy, zz when using split (linux)
    :param fastq_r1: name of the fastq file (read1)
    :param fname: the name of the file with no extensions, for naming the splits
    :param split_lines:  how many liens do you want each split file to be
    :return:  a list of the split names for the output of the split.smk rule
    """
    # All possible letter combinations
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
    
    name_list = total_names[:int(split_files)+1]  # slice only the names for the number of split files to be generated
    name_list2 = [split_fname.replace('_R1', '_R2') for split_fname in name_list]
    
    return name_list, name_list2


#############################################   Input and Variables  ##################################################

# get current username for computer
home_dir = str(Path.home())

# read in the config file to get pathways/samples into variables
configfile: f"{home_dir}/gbs/GBS_snakemake_pipeline/config/config.yaml"

read1 = config["sample_read1"]
read2 = read1.replace('_R1', '_R2')  # make the read2 name (only difference is R1 vs R2)
samplesheet = config["samplesheet"]
barcode_file = config["barcodefile"]

# get the pathway from the read1 fastq file, this is where we want to run the pipeline (set working dir)
path_name, base_name = os.path.split(read1)
file_handle = os.path.splitext(base_name)[0]
file_name = os.path.splitext(file_handle)[0]

# read the sample sheet and parse the data you need out of it into variables
index, samples, reflist, bams = parse_samplesheet(samplesheet)
ref_sample = zip(reflist, samples)

demuxlist, trimlist, bam_output, sort_output, unique_dirs = output_lists(ref_sample)
sam_ref = reference_samples(samples, reflist)
fastqc_output = [i.replace('sorted_bams', 'log/fastqc') for i in sort_output]

# how big do you want the chunk files to be.  default is 100M.  If memory is an issue, try smaller (like 10M).
split_number = 100000000  # must be a multiple of 4!!!!!

workdir: path_name

lib_id = file_name.split('_')[0]  # just get the library name without _R1 or _R2 for readgroup

#######################################  RULES ##################################################

rule all:
    input:
        expand("{sort_out}.bam", sort_out=sort_output),
        expand("{ref_dirs}/log/multiqc_report.html", ref_dirs=unique_dirs)

include: f"{home_dir}/gbs/GBS_snakemake_pipeline/workflow/rules/split.smk"
include: f"{home_dir}/gbs/GBS_snakemake_pipeline/workflow/rules/demultiplex.smk"
include: f"{home_dir}/gbs/GBS_snakemake_pipeline/workflow/rules/trimmomatic.smk"
include: f"{home_dir}/gbs/GBS_snakemake_pipeline/workflow/rules/bwa_map.smk"
include: f"{home_dir}/gbs/GBS_snakemake_pipeline/workflow/rules/novosort.smk"
include: f"{home_dir}/gbs/GBS_snakemake_pipeline/workflow/rules/multiqc.smk"
include: f"{home_dir}/gbs/GBS_snakemake_pipeline/workflow/rules/zipup.smk"
