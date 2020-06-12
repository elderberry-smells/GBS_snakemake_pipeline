rule zipup_demux:
    # we can zip up everything that remains (except for the bams which we need for VCF creation)
    input:
        multi = expand("{ref_dir}/log/multiqc_report.html", ref_dir=unique_dirs)
    output: expand("{ref_dir}/demultiplex.tar.gz", ref_dir=unique_dirs)
    threads: 1
    run:
        for i in unique_dirs:
            shell("tar -czvf {i}/demultiplex.tar.gz {i}/demultiplex/; rm -rf {i}/demultiplex/")


rule zipup_unmatched:
    input:
        # run after the demultiplex folders have been all zipped up
        demux_zip = expand("{ref_dir}/demultiplex.tar.gz", ref_dir=unique_dirs)
    output: "unmatched.tar.gz"
    threads: 1
    shell: "tar -czvf unmatched.tar.gz unmatched/ ; rm -rf unmatched/"
