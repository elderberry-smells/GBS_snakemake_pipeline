rule fastqc:
    input:
        # this input is to make sure that the sorted bams have been generated before finishing the stats
        fq1 = expand("sorted_bams/{sort_out}.bam", sort_out=sample_outputs),
        demux = expand("demultiplex/{sample}.{read}.fastq", sample=sample_outputs, read=[1,2])
    output:
        qc = expand("log/fastqc/{sample}.{read}_fastqc.html", sample=sample_outputs, read=[1,2])
    threads: 6
    shell:
         "fastqc --quiet -o log/fastqc/ -t {threads} {input.demux}"


rule multiqc:
    input:
        fastqc_files = expand("log/fastqc/{fastqc}.{read}_fastqc.html", fastqc=sample_outputs, read=[1, 2]),
        bams = expand("sorted_bams/{sort_out}.bam", sort_out=sample_outputs)
    output:
        "log/multiqc_report"
    shell:
          "multiqc log/ -o {output}"