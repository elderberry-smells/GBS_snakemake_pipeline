# VCF Analaysis Tools

| Call | Extended Call | Required | Value | Description |
| :--: | :----: | :-----: | :----------: | :--------------- |
| -vcf | --vcf | Yes | PATH | VCF file you want to filter.  Can be compressed with gzip or uncompressed |
| -out | --outfile | Yes | TEXT / PATH | Prefix for this filter, this will be the name of a new directory where the filytered VCF file will be saved (relative to path of VCF) |
| -qual | --min_quality | Optional | INT 0-999 | Set the QUAL minimum for the SNP site, remove SNP if score is below [qual] argument |
| -miss | --max_missing | Optional | FLOAT 0.0 - 1.0 | Remove SNP site if site doesn't have AT LEAST [max_missing] calls across all samples in site |
| -maf | --minor_allel_frequency | Optional | FLOAT 0.00 - 1.00 | "Remove SNP site if minor allele frequency is LOWER than [minor_allele_freq] |
| -mindp | --minimum_snp_depth | Optional | INT 1 - x | Remove SNP site when snp depth is BELOW [minimimum_snp_depth], over_rides homo/het depth if higher |
| -dphom | --min_depth_homozygous | Optional | INT 1 - x | Do not count homozygous SNP for sample if depth doesn't match or exceed [min_depth_homozygous] |
| -dphet | --min_depth_heterozygous | Optional | INT 1 - x | Do not count heterozygous SNP for sample if depth doesn't match or exceed [min_depth_heterozygous] |
| -removeindels | --removeindels | Optional | NO INPUT | Remove SNP site if REF or ALT is an indel. No value required after [-removeindels] |
| -csv | --csvout | Optional | NO INPUT | After filtering, remove headers and return a human readable CSV version of the VCF. No value required after [-csv] |
