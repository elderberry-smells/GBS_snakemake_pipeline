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


# read in the config file to get pathways/samples into variables
configfile: "~/gbs/GBS_snakemake_pipeline/config/config.yaml"
samplesheet = config["samplesheet"]
read1 = config["sample_read1"]

# get the pathway from the read1 fastq file, this is where we want to run the pipeline (set working dir)
path_name, base_name = os.path.split(read1)

# read the sample sheet and parse the data you need out of it into variables
index, samples= parse_samplesheet(samplesheet)

bams = bam_list(path_name, samples)  # create a list of bams to be run through mpileup at end stage

workdir: path_name

#######################################  RULES ##################################################

rule all:
    input:
        "calls/snps.raw.vcf.gz"  # final output is a single vcf file

include: "~/gbs/GBS_snakemake_pipeline/workflow/rules/demultiplex.smk"
include: "~/gbs/GBS_snakemake_pipeline/workflow/rules/trimmomatic.smk"
include: "~/gbs/GBS_snakemake_pipeline/workflow/rules/bwa_map.smk"
include: "~/gbs/GBS_snakemake_pipeline/workflow/rules/novosort.smk"
include: "~/gbs/GBS_snakemake_pipeline/workflow/rules/samtools_call.smk"
