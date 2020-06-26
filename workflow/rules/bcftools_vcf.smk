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
    output: 
        vcf = "filtered/{vcf_name}.recode.vcf",
        logfile = "filtered/{vcf_name}.log"
    wildcard_constraints:
        vcf_name= '[\w]*'
    threads: 1
    shell:
        "vcftools --vcf {input} "
        "--max-missing 0.7 "
        "--maf 0.02 "
        "--minQ 30 "
        "-recode --recode-INFO-all --out {params}"


rule summarize_filter:
    input:
         "filtered/{vcf_name}.log"
    params:
        vcf_dir = "chromosomes/filtered",
        scrip = f"{home_dir}/gbs/GBS_snakemake_pipeline/workflow/scripts/summary_vcf.py"
    output: "filtered/filter_summary.txt"
    wildcard_constraints:
        vcf_name= '[\w]*'
    shell:
        "python3 {params.scrip} -v {params.vcf_dir} -o {output} -s False"
