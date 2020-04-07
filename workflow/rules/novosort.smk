rule novosort:
    input:
        bam = "{ref_dir}/mapped_reads/{sample}.bam"
    output:
        sort_bam = "{ref_dir}/sorted_bams/{sample}.bam",
        index_bam = "{ref_dir}/sorted_bams/{sample}.bam.bai"
    log:  "{ref_dir}/log/novosort.log"
    threads: 16
    shell:
        "(novosort {input.bam} --threads {threads} --markDuplicates --index --output {output.sort_bam}) 2>> {log}"