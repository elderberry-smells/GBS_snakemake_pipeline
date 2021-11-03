rule novosort:
    input:
        bam = "mapped_reads/{sample}.bam"
    output:
        sort_bam = "sorted_bams/{sample}.bam",
        index_bam = "sorted_bams/{sample}.bam.bai"
    log:  "log/novosort.log"
    threads: 1
    shell:
        "(novosort {input.bam} --threads {threads} --markDuplicates --index --output {output.sort_bam}) 2>> {log}"
