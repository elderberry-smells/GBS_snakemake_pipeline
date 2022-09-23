rule zipup_demux:
    # we can zip up everything that remains (except for the bams which we need for VCF creation)
    input:
        mc = "log/multiqc_report",
        demux = f"{demux_parent}/"
        
    output: demux_output
    threads: 6
    shell:
         "tar cf - {input.demux} | pigz -9 -p {threads} > {output}; rm -rf {input.demux}; rm -rf log/fastqc"
