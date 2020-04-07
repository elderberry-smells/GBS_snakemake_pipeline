rule demultiplex:
    input:
        chunks1 = expand("chunks/{split_file}", split_file=split_R1),   # chunked read1 data
        chunks2 = expand("chunks/{split_file}", split_file=split_R2),   # chunked read2 data
        fq = read1,                                                     # path to fastq read1
        bar = barcode_file,                                             # path to barcode file
        samp = samplesheet                                              # path to samplesheet
    params:
        demux = f'{home_dir}/gbs/GBS_snakemake_pipeline/workflow/scripts/PE_fastq_demultiplex_AAFC.py'
    output:
        expand("{demux}.1.fastq", demux=demuxlist),
        expand("{demux}.2.fastq", demux=demuxlist)
    shell:
        "python3 {params.demux} -f {input.fq} -b {input.bar} -s {input.samp}"
