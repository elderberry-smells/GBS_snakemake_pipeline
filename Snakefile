def parse_samplesheet(sample_sheet):
    '''open samplesheet and grab relevant sample names.  Format of the sample sheet should be:
    sample #    index   sample_id   reference_path'''

    index = []
    samples = []
    ref_path = []
    with open(samplesheet, 'r') as samps:
        next(samps)  # skip the header line
        lines = samps.readlines()
        for line in lines:
            data = line.split()
            index.append(data[1])
            samples.append(data[2])
            ref_path.append(data[3])

    return index, samples, ref_path


def bam_list(sample_list):
    "make a list of bam files to run through the samtools mpileup"
    with open("data/sample3/bams.list", "w") as bam:
        for i in sample_list:
            line = f'data/sample3/sorted_reads/{i}.sorted.bam\n'
            bam.write(line)

####    initialize all the data into variables from the samplesheet
samplesheet = 'data/sample3/sample3_samplesheet.txt'
index, samples, ref_path = parse_samplesheet(samplesheet)  # returns a list of each important column in samplesheet
bam_list(samples)  # create a list of bams to be run through mpileup

#workdir:"PycharmProjects/gbs_pipeline/data/sample3"

rule all:
    input:
        # final output is a single vcf file
        "data/sample3/calls/snps.raw.vcf.gz"
       #'/calls/snps.raw.vcf.gz'


rule demultiplex:
    input:
        fq = "data/sample3/sample3_R1.fastq",
        bar = "/home/bioinf/PycharmProjects/gbs_pipeline/barcodes/barcodes_192.txt",
        samp = samplesheet
    output:
        expand("data/sample3/demultiplex/{sample}.1.fastq", sample=samples),
        expand("data/sample3/demultiplex/{sample}.2.fastq", sample=samples)
    conda:
    shell:
        "python3 scripts/fastq_demultiplex3.py -f {input.fq} -b {input.bar} -s {input.samp}"


rule trimmomatic:
    input:
        f1 = "data/sample3/demultiplex/{sample}.1.fastq",
        f2 = "data/sample3/demultiplex/{sample}.2.fastq"
    output:
        p1 = "data/sample3/trimmomatic/{sample}.1.paired.fastq",
        u1 = "data/sample3/trimmomatic/{sample}.1.unpaired.fastq",
        p2 = "data/sample3/trimmomatic/{sample}.2.paired.fastq",
        u2 = "data/sample3/trimmomatic/{sample}.2.unpaired.fastq"
    log: "data/sample3/log/trimmomatic.output.log"
    threads: 4
    shell:  "(trimmomatic PE -threads {threads} -phred33 {input.f1} {input.f2} {output.p1} {output.u1} {output.p2} {output.u2} "
            "ILLUMINACLIP:scripts/trimmomatic-0.39-1/TruSeq2-PE.fa:2:40:15 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:55) 2>> {log}"


rule bwa_map:
    input:
        ref = "reference/barley/160404_barley_pseudomolecules_masked.fasta",
        r1 = "data/sample3/trimmomatic/{sample}.1.paired.fastq",
        r2 = "data/sample3/trimmomatic/{sample}.2.paired.fastq"
    output:
        bam = "data/sample3/mapped_reads/{sample}.bam"
        #bam = temp("data/sample3/mapped_reads/{sample}.bam")
    log: "data/sample3/log/bwa_map.log"
    threads: 4
    shell: "(bwa mem {input.ref} {input.r1} {input.r2} | samtools view -Shbu > {output.bam}) 2>> {log}"


rule novosort:
    input:
        bam = "data/sample3/mapped_reads/{sample}.bam"
    output:
        sort_bam = "data/sample3/sorted_reads/{sample}.sorted.bam",
        index_bam = "data/sample3/sorted_reads/{sample}.sorted.bam.bai"
    log:
        "data/sample3/log/novosort.log"
    threads: 4
    shell:
        "novosort {input.bam} --threads {threads} --index --output {output.sort_bam} 2>> {log}"


rule mpileup_call:
    input:
        expand("data/sample3/sorted_reads/{sample}.sorted.bam", sample=samples),
        ref = "reference/barley/160404_barley_pseudomolecules_masked.fasta",
        bammies = temp("data/sample3/bams.list")
    output:
        "data/sample3/calls/snps.raw.vcf.gz"
    threads: 4
    shell:
        "samtools mpileup -u -t AD,DP -f {input.ref} -b {input.bammies} | bcftools call -mv -Oz > {output}"


#  snakemake --cluster "qsub -pe threaded {threads}" --jobs 100