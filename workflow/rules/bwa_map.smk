rule bwa_map:
    input:
        ref = ref_genome,
        r1 = "trimmomatic/{sample}.1.paired.fastq",
        r2 = "trimmomatic/{sample}.2.paired.fastq"
    output:
        bam = temp("mapped_reads/{sample}.bam")
    log: "log/bwa_map.log"
    threads: 16
    shell: "(bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | samtools view -Shbu > {output.bam}) 2>> {log}"
