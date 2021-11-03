rule fastqc:
    input:
        # this input is to make sure that the sorted bams have been generated before finishing the stats
        "demultiplex/{sample}.{read}.fastq"
    output:
        "log/fastqc/{sample}.{read}_fastqc.html"
    threads: 1
    shell:
         "fastqc -o log/fastqc/ -t {threads} {input}"
