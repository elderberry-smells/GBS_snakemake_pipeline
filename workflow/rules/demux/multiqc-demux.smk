rule multiqc:
    input:
        fastqc_files = expand("log/fastqc/{fastqc}.{read}_fastqc.html", fastqc=fastqc_output, read=[1, 2]),
        bams = expand("{sort_out}.bam", sort_out=sort_output)
    output:
        "log/multiqc_report.html"
    run:
        for i in unique_dirs:
            shell("multiqc {i}/log -o {i}/log/")
