rule zipup_demux:
    # we can zip up everything that remains (except for the bams, those will need to be set until VCF made)
    input: 
        fastqc_files = expand("{reference_folder}/{fastqc}.{read}_fastqc.html", reference_folder=unique_dirs, fastqc=fastqc_output, read=[1, 2]),
        demux_folder = "{reference_folder}/demultiplex/"  # get all the demux folders from pipeline
    output: "{reference_folder}/demultiplex.tar.gz"
    threads: 1
    shell: "tar -czvf {output} {input.demux_folder}; rm -rf {input.demux_folder}"


rule zipup_unmatched:
    input: 
        demux_zip = expand("{reference_folder}/demultiplex.tar.gz", reference_folder=unique_dirs),
        unmatched_folder = "unmatched/"
    output: 'unmatched.tar.gz'
    threads: 1
    shell:  "tar -czvf {output} {input.unmatched_folder}; rm -rf {input.unmatched_folder}"

