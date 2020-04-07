split_R1, split_R2 = generate_split_list(read1, file_name, split_number)  # get a list of split output
chunk_prefix1 = f"chunks/{file_name}_"
chunk_prefix2 = chunk_prefix1.replace('_R1', '_R2')

rule split:
    input:
        r1 = read1,
        r2 = read2
    params:
        p1 = chunk_prefix1,
        p2 = chunk_prefix2,
        l = split_number
    output:
        expand("chunks/{split_file}", split_file=split_R1),
        expand("chunks/{split_file}", split_file=split_R2)
    shell:
        "zcat {input.r1} | split -l {params.l} --additional-suffix=.fq - {params.p1}; "
        "zcat {input.r2} | split -l {params.l} --additional-suffix=.fq - {params.p2}"
        
