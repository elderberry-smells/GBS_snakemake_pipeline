rule bwa_map:
    input:
        r1 = "trimmomatic/{sample}.1.paired.fastq",
        r2 = "trimmomatic/{sample}.2.paired.fastq"
    output:
        bam = "mapped_reads/{sample}.bam"
    params:
        readgroup = lambda wildcards: f"@RG\\tID:{lib_id}\\tSM:{wildcards.sample}\\tPL:ILLUMINA",
        ref = lambda wildcards: sam_ref[wildcards.sample]
    log: "log/bwa_map.log"
    threads: 1
    shell:
        "(bwa mem -R '{params.readgroup}' -t {threads} {params.ref} {input.r1} {input.r2} | samtools view -Shbu - > {output.bam}) 2>> {log}"
