# VCF Analaysis Tools

## **1 . VCF Filter Program**

### **Description**
This is a filtering script for VCF files (VCF format version 4.2) written in python (3.6).  This script allows you to set parameters and filter out SNP sites that match or exceed parameters indicated. 

Filter VCF on missing, quality, depth (general depth, depth for homozygous calls, depth for heterozygous calls), indels (remove or retain) and minor allele frequency.

### **Inputs and Outputs for Program**

- VCF file can be compressed, or uncompressed.  

- Output is designated as a directory, and the filtered VCF will save there, as well as some stats on the VCF filter.

- CSV output  
    - add -csv to arguments to save the output as a CSV format for transfer to windows or other machine
    - best suited to run this program once you are satisfied with final filtering of VCF, as this program can not filter a CSV file. This also automatically removes indels.


### **How to use the filter progam:**
Navigate to the directory of the VCF file and call the program with thte following:

```
python3 ~/gbs/GBS_snakemake_pipeline/workflow/scripts/AAFC_filter_vcf.py 
```
Following the command, you can add your arguments for filtering, as seen in options below.


### **Arguments for filter program:**
| Call | Extended Call | Required | Value | Description |
| :--: | :----: | :-----: | :----------: | :--------------- |
| -vcf | --vcf | Yes | PATH | VCF file you want to filter.  Can be compressed with gzip or uncompressed |
| -out | --outfile | Yes | TEXT / PATH | Prefix for this filter, this will be the name of a new directory where the filtered VCF file will be saved (relative to path of CWD) |
| -qual | --min_quality | Optional | INT 0-999 | Set the QUAL minimum for the SNP site, remove SNP if score is below [qual] argument |
| -miss | --max_missing | Optional | FLOAT 0.0 - 1.0 | Remove SNP site if site doesn't have AT LEAST [max_missing] calls across all samples in site. Written as % of samples to have a call ie. 0.7 == 70% of samples with call) |
| -maf | --minor_allel_frequency | Optional | FLOAT 0.00 - 1.00 | "Remove SNP site if minor allele frequency is LOWER than [minor_allele_freq] |
| -mindp | --minimum_snp_depth | Optional | INT 1 - x | Remove SNP site when snp depth is BELOW [minimimum_snp_depth], over_rides homo/het depth if higher |
| -dphom | --min_depth_homozygous | Optional | INT 1 - x | Do not count homozygous SNP for sample if depth doesn't match or exceed [min_depth_homozygous] |
| -dphet | --min_depth_heterozygous | Optional | INT 1 - x | Do not count heterozygous SNP for sample if depth doesn't match or exceed [min_depth_heterozygous] |
| -removeindels | --removeindels | Optional | NO INPUT | Remove SNP site if REF or ALT is an indel. No value required after [-removeindels] |
| -csv | --csvout | Optional | NO INPUT | After filtering, remove headers and return a human readable CSV version of the VCF. No value required after [-csv] |

### **Stats output for the filter program**
The program will have 2 outputs when run, both of which will be located in the -out directory in your arguments.

- Filtered VCF - This file will be named the same as your input VCF file, but with a .filtered.vcf suffix
- Stas for filter - This text file will show stats for the VCF file after filtering, which includes how many sites remaining, how many sites were removed due to filter aruments, and the general stats per chromosome in your VCF.  

An example of a stats output is shown below, which came from a command line call as follows:

```
$ python3 ~/gbs/GBS_snakemake_pipeline/workflow/scripts/AAFC_filter_vcf.py -vcf example_vcf.vcf.gz -out filter/standard_filter -qual 30 -miss 0.7 - maf 0.05 -dphom 1 -dphom 2 -removeindels
```
filter/standard_filter/example_vcf.stats
```
AAFC VCF FILTER
Author: Brian James (brian.james4@canada.ca)
Created: SEPT 2020

Parameters as interpreted:
	--vcf example_vcf.vcf.gz
	--out filtered/standard_filter/example_vcf.filtered.vcf
	--qual 30
	--miss 0.7
	--maf 0.05
	--dphom 1
	--dphet 2
	--removeindels
	

VCF FILTER STATS:
After filtering, kept 19639 out of a possible 9416529 sites
Total Samples in VCF: 1477

SNP sites dropped due to the following parameters:
	QUAL:	4486604
	INDEL:	193849
	DP:	46480
	MISS:	4494943
	MAF:	175014

Genome Distribution of Remaining SNP Sites and Stats:

CHROM	SITES	#SNPs	AVG.DP	%HET	AVG.MAF
chr1H	2406	2871473	1.96	18.47	0.285
chr2H	3089	3681786	1.93	19.30	0.303
chr3H	2875	3407651	1.84	17.65	0.283
chr4H	1843	2212995	1.88	17.64	0.261
chr5H	2902	3432777	1.83	17.58	0.278
chr6H	2459	2925840	1.87	18.91	0.274
chr7H	3675	4427481	2.00	19.34	0.282
chrUn	390	473668	2.05	25.12	0.294            
```

## **2. Windowed Variance**

### **Description**
Analyze a VCF file (VCF format version 4.2) for # of SNPs across genome, in a windowed breakdown of each chromosome. The program will output a colored graph that will plot the number of SNPs per window across each chromosome.

The VCF file can be uncompressed or compressed (gzip) and **MUST** have headers for this program to work!

The headers contain the information about each chromosomes (names and lengths) required for the graphing and breakdown. 

### **How to use the program**
Simply direct the program to a VCF file and the window size you would like to analyze.  Window is just a number of positions along the genome to sum the # of SNPs in.  

The output will be saved into the same directory as the VCF file, and be named vcf_name.coverage.png

example:
```
$ python3 ~/gbs/GBS_snakemake_pipeline/workflow/scripts/windowed_variance.py -v example_vcf.vcf -w 1000
```

### Graphical Output of Program
The following output is from an output of over ~1500 barley lines in a VCF file. The genome consists of 7 chromosomes and Scaffolds (designated CHRUN)
![windowed_variance](https://github.com/elderberry-smells/GBS_snakemake_pipeline/blob/master/workflow/resources/images/example_VCF_coverage.png?raw=true)

