rule fastqc:
    input:
        sbams = expand("{sort_out}.bam", sort_out=sort_output)
    output:
        expand("{fastqc}.{read}_fastqc.html", fastqc=fastqc_output, read=[1, 2])
    threads: 16
    run:
        for i in unique_dirs:
            shell("fastqc -o {i}/log/fastqc/ -t {threads} {i}/demultiplex/*.fastq")


rule multiqc:
    input:
        expand("{fastqc}.{read}_fastqc.html", fastqc=fastqc_output, read=[1, 2])
    output:
        expand("{ref_dir}/log/multiqc_report.html", ref_dir=unique_dirs)
    run:
        for i in unique_dirs:
            shell("multiqc {i}/log -o {i}/log/")
