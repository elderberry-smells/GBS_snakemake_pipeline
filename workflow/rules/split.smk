split_R1, split_R2 = generate_split_list(read1, file_name, split_number)  # get a list of split output
chunk_prefix1 = f"chunks/{file_name}_"

# formatting for if read number is lowercase or uppercase
if '_r1' in chunk_prefix1:
    chunk_prefix2 = chunk_prefix1.replace('_r1', '_r2')
else:
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
        temp(expand("chunks/{split_file}", split_file=split_R1)),
        temp(expand("chunks/{split_file}", split_file=split_R2))
    shell:
        "zcat {input.r1} | split -l {params.l} --additional-suffix=.fq - {params.p1}; "
        "zcat {input.r2} | split -l {params.l} --additional-suffix=.fq - {params.p2}"
        
