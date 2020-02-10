rule novosort:
    input:
        bam = "mapped_reads/{sample}.bam"
    output:
        sort_bam = "sorted_reads/{sample}.sorted.bam",
        index_bam = "sorted_reads/{sample}.sorted.bam.bai"
    log:  "log/novosort.log"
    threads: 4
    shell:
        "novosort {input.bam} --threads {threads} --index --output {output.sort_bam} 2>> {log}"