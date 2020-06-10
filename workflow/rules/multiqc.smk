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
        # make the zip up a requirement so everything gets compressed before finishing the stats
        zip_demux = expand("{reference_folder}/demultiplex.tar.gz", reference_folder=unique_dirs),
        zip_unmatched = "unmatched.tar.gz",
        expand("{fastqc}.{read}_fastqc.html", fastqc=fastqc_output, read=[1, 2])
    output:
        expand("{ref_dir}/log/multiqc_report.html", ref_dir=unique_dirs)
    run:
        for i in unique_dirs:
            shell("multiqc {i}/log -o {i}/log/")
