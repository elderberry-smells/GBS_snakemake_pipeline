rule fastqc:
    input:
        # this input is to make sure that the sorted bams have been generated before finishing the stats
        fq1 = expand("{sort_out}.bam", sort_out=sort_output)
    output:
        qc = "{sample}.{read}_fastqc.html"
    threads: 1
    shell:
         "fastqc -o log/fastqc/ -t {threads} demultiplex/*.fastq"


rule multiqc:
    input:
        fastqc_files = expand("{fastqc}.{read}_fastqc.html", fastqc=fastqc_output, read=[1, 2]),
        bams = expand("{sort_out}.bam", sort_out=sort_output)
    output:
        "log/multiqc_report"
    shell:
          "multiqc log/ -o {output}"