rule demultiplex:
    input:
        fq = config["sample_read1"],
        bar = config["barcodefile"],
        samp = samplesheet
    output:
        expand("demultiplex/{sample}.1.fastq", sample=samples),
        expand("demultiplex/{sample}.2.fastq", sample=samples)
    shell:
        "python3 workflow/scripts/PE_fastq_demultiplex.py -f {input.fq} -b {input.bar} -s {input.samp}"