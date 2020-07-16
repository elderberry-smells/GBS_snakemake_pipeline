# seperate VCFs per chromosome will be generated
rule bcf_split:
    input:
        bammies = 'bams.list',
        ref = reference
    output:
        vcf = "chromosomes/{vcf_name}.{chrom}.vcf"
    wildcard_constraints:
        vcf_name = '[\w]*'
    threads: 1
    run:
        if wildcards.chrom != 'scaffolds':  # generate a VCF file for each chromosome
            shell("bcftools mpileup -a AD,DP -r {wildcards.chrom} -f {input.ref} -b {input.bammies} | bcftools call -mv -Ov > {output.vcf}")

        else:  # generate one VCF with all scaffolds
            shell("bcftools mpileup -a AD,DP -R scaffolds.txt -f {input.ref} -b {input.bammies} | bcftools call -mv -Ov > {output.vcf}")


rule merge_vcf:
    input:  expand("chromosomes/{vcf_name}.{chrom}.vcf", vcf_name=p_dir, chrom=chr_list)
    output:  "chromosomes/{vcf_name}.vcf"
    wildcard_constraints:
        vcf_name = '[\w]*'
    threads: 1
    shell:  "bcftools concat chromosomes/*.vcf > {output}"
