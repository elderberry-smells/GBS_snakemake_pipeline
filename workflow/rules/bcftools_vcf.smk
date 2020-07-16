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
