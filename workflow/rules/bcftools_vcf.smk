# the entire VCF will be produced and filtered
rule bcf_split:
    input:
        bammies = 'bams.list',
        ref = reference
    output:  "{vcf_name}.vcf"
    wildcard_constraints:
        vcf_name = '[\w]*'
    threads: 1
    shell:
          "bcftools mpileup -a AD,DP -f {input.ref} -b {input.bammies} | bcftools call -mv -Ov > {output}"


# add in filtering of the VCF file for specific quality parameters
rule filter_vcf:
    input: "{vcf_name}.vcf"
    params: "filtered/{vcf_name}"
    output: "filtered/{vcf_name}.recode.vcf"
    wildcard_constraints:
        vcf_name= '[\w]*'
    threads: 1
    shell:
        "vcftools --vcf {input} "
        "--max-missing 0.7 "
        "--maf 0.02 "
        "--minQ 30 "
        "-recode --recode-INFO-all --out {params}"

