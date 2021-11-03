rule trimmomatic:
    input:
        f1 = "demultiplex/{sample}.1.fastq",
        f2 = "demultiplex/{sample}.2.fastq"
    
    output:
        p1 = "trimmomatic/{sample}.1.paired.fastq",
        u1 = "trimmomatic/{sample}.1.unpaired.fastq",
        p2 = "trimmomatic/{sample}.2.paired.fastq",
        u2 = "trimmomatic/{sample}.2.unpaired.fastq"
    params:
        trim = f"{home_dir}/anaconda3/envs/gbs/share/trimmomatic-0.39-1/adapters/TruSeq2-PE.fa"
    log:
        "log/trimmomatic.output.log"
    threads: 1
    shell:  "(trimmomatic PE -threads {threads} -phred33 {input.f1} {input.f2} " \
            "{output.p1} {output.u1} {output.p2} {output.u2} " \
            "ILLUMINACLIP:{params.trim}:2:40:15 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:55) 2>> {log}"
