rule trimmomatic:
    input:
        f1 = trim_input1,
        f2 = trim_input2   
    output:
        p1 = temp("trimmomatic/{sample}.1.paired.fastq"),
        u1 = temp("trimmomatic/{sample}.1.unpaired.fastq"),
        p2 = temp("trimmomatic/{sample}.2.paired.fastq"),
        u2 = temp("trimmomatic/{sample}.2.unpaired.fastq")
    params:
        trim = f"{home_dir}/miniconda3/envs/gbs/share/trimmomatic-0.39-1/adapters/TruSeq2-PE.fa"
    log:
        "log/trimmomatic.output.log"
    threads: 1
    shell:  "(trimmomatic PE -threads {threads} -phred33 {input.f1} {input.f2} " \
            "{output.p1} {output.u1} {output.p2} {output.u2} " \
            "ILLUMINACLIP:{params.trim}:2:40:15 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:55) 2>> {log}"
