rule demultiplex:
    input:
        fq = read1,
        bar = barcode_file,
        samp = samplesheet
    output:
        expand("demultiplex/{sample}.1.fastq", sample=samples),
        expand("demultiplex/{sample}.2.fastq", sample=samples)
    shell:
        "python3 /home/AAFC-AAC/your_user_name/gbs/GBS_snakemake_pipeline/workflow/scripts/PE_fastq_demultiplex.py -f {input.fq} -b {input.bar} -s {input.samp}"
