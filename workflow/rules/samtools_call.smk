rule samtools_call:
    input:
        expand("sorted_reads/{sample}.sorted.bam", sample=samples),
        ref = config["reference_file"],
        bammies = bams
    output:
        "calls/snps.raw.vcf.gz"
    threads: 16
    shell:
        "samtools mpileup -u -t AD,DP -f {input.ref} -b {input.bammies} | bcftools call -mv -Oz > {output}"