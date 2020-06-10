rule zipup_demux:
    # we can zip up everything that remains (except for the bams, those will need to be set until VCF made)
    input: "{reference_folder}/demultiplex/"  # get all the demux folders from pipeline
    output: "{reference_folder}/demultiplex.tar.gz"
    threads: 1
    shell: "tar -czvf {output} {input}; rm -rf {input}"


rule zipup_unmatched:
    input: 'unmatched/'
    output: 'unmatched.tar.gz'
    threads: 1
    shell:  "tar -czvf {output} {input}; rm -rf {input}"

