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
        if wildcards.chrom != 'scaffolds':
            shell("bcftools mpileup -a AD,DP -r {wildcards.chrom} -f {input.ref} -b {input.bammies} "
                  "| bcftools call -mv -Ov > {output.vcf}")
        else:
            shell("bcftools mpileup -a AD,DP -R scaffolds.txt -f {input.ref} -b {input.bammies} "
                  "| bcftools call -mv -Ov > {output.vcf}")


# filtering of the VCF file for specific quality parameters
rule filter_vcf:
    input: "chromosomes/{vcf_name}.{chrom}.vcf"
    params:
        vcf = "chromosomes/filtered/{vcf_name}.{chrom}"
    output:
        filtered = "chromosomes/filtered/{vcf_name}.{chrom}.recode.vcf",
        log = "chromosomes/filtered/{vcf_name}.{chrom}.vcf.log"
    wildcard_constraints:
        vcf_name= '[\w]*'
    threads: 1
    shell:
        "vcftools --vcf {input} "
        "--max-missing 0.7 "
        "--minQ 30 "
        "--recode --recode-INFO-all --out {params.vcf}"
