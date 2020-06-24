rule zipup_demux:
    # we can zip up everything that remains (except for the bams which we need for VCF creation)
    input:
        multi = expand("{ref_dir}/log/multiqc_report.html", ref_dir=unique_dirs)
    output: expand("{ref_dir}/{ref_dir}_demux.tar.gz", ref_dir=unique_dirs)
    threads: 16
    run:
        for i in unique_dirs:
            shell("tar cf - {i}/demultiplex/ | pigz -9 -p {threads} > {i}/{i}_demux.tar.gz; rm -rf {i}/demultiplex/")


rule zipup_unmatched:
    input:
        # run after the demultiplex folders have been all zipped up
        demux_zip = expand("{ref_dir}/{ref_dir}_demux.tar.gz", ref_dir=unique_dirs)
    output: "unmatched.tar.gz"
    threads: 16
    shell: "tar cf - unmatched/ | pigz -9 -p {threads} > unmatched.tar.gz; rm -rf unmatched"
