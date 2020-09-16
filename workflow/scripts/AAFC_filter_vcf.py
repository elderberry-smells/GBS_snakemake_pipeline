#!/usr/bin/env python3
from pathlib import Path
import csv
import re
import argparse
import gzip


def get_arguments():

    parser = argparse.ArgumentParser(description='Get a summary of vcf file depth')
    parser.add_argument('-vcf', '--vcf', required=True, help="VCF file you want to filter.  Can be compressed with gzip or uncompressed")
    parser.add_argument('-out', '--outfile', required=True, help="Name and path for filtered VCF file")
    parser.add_argument('-csv', '--csvout', action='store_true', help='After filtering, remove headers and return a human readable CSV version of the VCF')
    parser.add_argument('-qual', '--min_quality', default=False, help='Set the QUAL minimum for the SNP site, remove SNP if score is below [qual] argument')
    parser.add_argument('-miss', '--max_missing', default=False, help="Remove SNP site if site doesn't have AT LEAST [max_missing] calls across pop.")
    parser.add_argument('-maf', '--minor_allele_freq', default=False, help="Remove SNP site if minor allele frequency is lower than [minor_allele_freq]")
    parser.add_argument('-removeindels', '--removeindels', action='store_true', help="Remove SNP site if REF or ALT is an indel")
    parser.add_argument('-dphom', '--min_depth_homozygous', default=False, help="Do not count homozygous SNP for sample if depth doesn't match or exceed [min_depth_homozygous]")
    parser.add_argument('-dphet', '--min_depth_heterozygous', default=False, help="Do not count heterozygous SNP for sample if depth doesn't match or exceed [min_depth_homozygous]")
    parser.add_argument('-mindp', '--minimimum_snp_depth', default=False, help="Remove SNP site when snp depth is BELOW [minimimum_snp_depth], over_rides homo/het depth if higher")
    args = parser.parse_args()
    return args


def generate_outfile(infile_path, outfile_params, csv_output):
    # take the arguments for infile/outfile and generate the directory and paths for the filter out and stats out
    # figure out here if the file or directory already exists, then inside the directory add the (1) or (2) or whatever copy this would be.
        
    # take the out argument and create a directory with that prefix   
    outdir = Path(outfile_params)
    outdir.mkdir(parents=True, exist_ok=True)
    
    if infile_path.suffix == '.gz':
        vcf_name = infile_path.stem
        vcf = Path(vcf_name).stem
    else:
        vcf = Path(infile_path).stem
    
    if csv_output:
        vcf_out = Path(outdir / f'{vcf}.filter.csv')
    else:
        vcf_out = Path(outdir / f'{vcf}.filter.vcf')
    stats_out = Path(outdir / f'{vcf}.stats')
    
    return vcf_out, stats_out


def filter_vcf(infile, filter_params, outfiles):
  
    summarized_data = {}
    sample_list = []
    chromosomes = []
    population = 0
    total_snps = 0
    post_filter = 0
    
    why = {'QUAL': 0, 'INDEL': 0, 'DP': 0, 'MISS': 0, 'MAF': 0}
    
    for line in infile:
        
        try:
            data = line.decode("utf-8").strip()  # decode because VCF was compressed and is in binary mode
        except:
            data = line.strip()  # just strip the line, no compresion in file
        
        if data.startswith('##'):
            if 'contig' in data.lower():  # isolates the header lines detailing the chromosomes and scaffold information
                
                if 'scaffold' not in data.lower():
                    chrom_name = re.findall(r'ID=(.*),', data)[0]
                    
                    # add chromosome name to summarized_data and make it an empty dict to add samples to
                    chromosomes.append(chrom_name)
                    # write the header line to the ouput file (if not csv)
                    if filter_params['CSV']:
                        continue
                    else:
                        try:
                            outfiles.writelines(line)
                        except:
                            outfiles.writelines(line.decode("utf-8"))
                            
                else:
                    if 'scaffold' in data.lower():
                        # add a general "scaffold" to the dictionary for stats purposes
                        # if 'scaffold' in summarized_data.keys():
                        if 'scaffold' in chromosomes:
                            continue
                        else:
                            # summarized_data['scaffold'] = {}
                            chromosomes.append('Scaf')
                    
                        # write the header line to the ouput file (if not csv)
                        if filter_params['CSV']:
                            continue
                        else:
                            try:
                                outfiles.writelines(line)
                            except:
                                outfiles.writelines(line.decode("utf-8"))
                        
            else:  # remainder of header info that isn't contig details
                # write the header line to the ouput file (if not csv) 
                if filter_params['CSV']:
                    continue
                else:
                    try:
                        outfiles.writelines(line)
                    except:
                        outfiles.writelines(line.decode("utf-8"))  
        
        elif data.startswith("#CHROM"):
            # make a sample list and add up the population numbers
            for samp in data.split()[9:]:
                sample_list.append(samp)
                population += 1  # add the sample # to the population list
                
            # add the sample_list items to the dictionary
            for chrom in chromosomes:
                summarized_data[chrom] = {'sites': 0, 'snps':0, 'dp': 0, 'hets': 0, 'maf': [0, 0]}
                
            # write the line to the outfile
            if filter_params['CSV']:
                skip_col = ['INFO', 'FORMAT']
                header = [i for i in line.split() if i not in skip_col]
                outfiles.writerow(header)
            else:
                try:
                    outfiles.writelines(line)
                except:
                    outfiles.writelines(line.decode("utf-8"))
                                    
        else:
            
            total_snps += 1  # add to the total SNPs in VCF (pre-filter)
            snp = data.split()
            
            # initialize the data to write to the outfile with info columns
            if filter_params['CSV']:
                # skip info line and format line
                line_data = line_data = [snp[0], snp[1], snp[2], snp[3], snp[4], snp[5], snp[6]] 
            else:
                line_data = [snp[0], snp[1], snp[2], snp[3], snp[4], snp[5], snp[6], snp[7], snp[8]]  
            
            # vcf format for each line/site in file for easier reference
            chrom = snp[0]
            ref = snp[3]
            alt = snp[4]
            qual = float(snp[5])
            info = snp[7]
                        
            # screen site based on INFO before moving into specific format info for each sample (if info set)
            if filter_params['QUAL'] and qual < int(filter_params['QUAL']):
                why['QUAL'] += 1
                continue
            
            elif filter_params['INDEL'] and 'INDEL' in info:
                why['INDEL'] += 1
                continue
            
            else:
                # made it through the first 2 INFO based parameters, now filter based on sample info
            
                # parameters remaining to filter per sample/population: 'MISS', 'MAF', 'MIN_DP', 'DP_HOMO', 'DP_HET' 

                snp_site = {'ref': 0, 'alt':0, 'miss': 0, 'het':0, 'dp': [0, 0]}
                
                # determine if the depth is causing the dropped lines by adding to this dict
                dp_drops = {"DP": 0, "DP_HOMO": 0, "DP_HET": 0}
                                    
                for site in snp[9:]:  # just the sample information in a loop, get the values and if all match at end, write the line to the outfile
                    
                    gt, pl, dp, ad = site.split(':')
                    
                    # look into each sample, determine if the general parameters are met (general depth of sample)
                    if filter_params['MIN_DP'] and int(dp) < int(filter_params['MIN_DP']):
                        if gt != './.': # only add it to being a depth issue if called data
                            dp_drops['DP'] += 1
                        site = f"./.:{pl}:{dp}:{ad}"  # didn't pass the test, make it no data and append it to snp_site for tally
                        snp_site['miss'] += 1
                        if filter_params['CSV']:
                            line_data.append('NA')
                        else:
                            line_data.append(site)

                    else:
                        if gt == './.':  # missing data point
                            snp_site['miss'] += 1
                            # add the info to the line to write(even though it may not be written in the end after filtering for MISS)
                            if filter_params['CSV']:
                                line_data.append('NA')
                            else:
                                line_data.append(site)
                            
                        elif gt == '0/0':  # reference base called
                            if filter_params['DP_HOMO'] and int(dp) < int(filter_params['DP_HOMO']):
                                dp_drops['DP_HOMO'] += 1
                                site = site.replace('0/0:', './.:')  # didn't pass the test, make it no data and append it to snp_site for tally
                                snp_site['miss'] += 1
                                if filter_params['CSV']:
                                    line_data.append('NA')
                                else:
                                    line_data.append(site)
                                snp_site['miss'] += 1
                                
                            else:  # passed the homozygous filter for ref
                                snp_site['ref'] += 2
                                snp_site['dp'][0] += 1  # add snp # tally
                                snp_site['dp'][1] += int(dp)  # add depth to make avg. depth at end
                                if filter_params['CSV']:
                                    line_data.append(f'{ref}/{ref}')
                                else:
                                    line_data.append(site)                        
                            
                        elif gt == '1/1':  # reference base called
                            if filter_params['DP_HOMO'] and int(dp) < int(filter_params['DP_HOMO']):
                                dp_drops['DP_HOMO'] += 1
                                site = site.replace('1/1:', './.:')  # didn't pass the test, make it no data and append it to snp_site for tally
                                if filter_params['CSV']:
                                    line_data.append('NA')
                                else:
                                    line_data.append(site)
                                snp_site['miss'] += 1
                                               
                            else:  # passed the homozygous filter for ref
                                snp_site['alt'] += 2
                                snp_site['dp'][0] += 1  # add snp # tally
                                snp_site['dp'][1] += int(dp)  # add depth to make avg. depth at end
                                if filter_params['CSV']:
                                    line_data.append(f'{alt}/{alt}')
                                else:
                                    line_data.append(site)
                                                    
                        else:  # its a heterozygous SNP
                            if filter_params['DP_HET'] and int(dp) < int(filter_params['DP_HET']):
                                dp_drops['DP_HET'] += 1
                                site = site.replace('0/1:', './.:')  # didn't pass the test, make it no data and append it to snp_site for tally
                                
                                if filter_params['CSV']:
                                    line_data.append('NA')
                                else:
                                    line_data.append(site)
                                    snp_site['miss'] += 1
                                                           
                            else:  # passed the heterozygous filter for ref
                                snp_site['ref'] += 1
                                snp_site['alt'] += 1
                                snp_site['het'] += 1
                                snp_site['dp'][0] += 1  # add snp # tally
                                snp_site['dp'][1] += int(dp)  # add depth to make avg. depth at end
                                
                                if filter_params['CSV']:
                                    line_data.append(f'{ref}/{alt}')
                                else:
                                    line_data.append(site)                                                                  
                    
                # determine if the numbers for missing and MAF are met after looping through the samples data.  If yes, write the line, if not, skip.
                missing = 1 - (float(snp_site['miss']/population))
                try: # might be mostly missing data in a SNP and cause a zero division error.  This will be eliminated by missing param anyways.
                    allele_counts = [snp_site['ref'], snp_site['alt']]
                    major = float(max(allele_counts))
                    minor = float(min(allele_counts))
                    maf = 1 - (major/(major + minor))
                except:
                    maf = 0
                
                #filter on missing parameter (anything that is UNDER the parameter fails)
                if filter_params['MISS'] and float(missing) < float(filter_params['MISS']): 
                    
                    # determine if the the number of lines would still pass missing if these dp drops didn't happen:
                    if filter_params['MIN_DP'] or filter_params['DP_HOMO'] or filter_params['DP_HET']:   
                        
                        dp_total_drops = (dp_drops['DP'] + dp_drops['DP_HOMO'] + dp_drops['DP_HET'])         
                        no_dp_drops = int(snp_site['miss']) - int(dp_total_drops)
                        
                        if (1- float(no_dp_drops/population)) > float(filter_params['MISS']): 
                            # The parameter would have passed when you disregard dropped DP SNPs
                            why['DP'] += 1
                        else:
                            # the parameter failed regardless of SNP depth, and was strictly because of missing SNP data
                            why['MISS'] += 1
                    else:
                        why['MISS'] += 1 
                
                # filter on MAF parameter (anything that is UNDER the parameter fails)
                elif filter_params['MAF'] and float(maf) < float(filter_params['MAF']):
                    why['MAF'] += 1
                                                              
                else:
                    # the SNP site has passed the parameters to be accepted (add site to post filter) and write that sucker out
                    post_filter += 1
                    
                    if filter_params['CSV']:
                        outfiles.writerow(line_data)
                    else:
                        outfiles.write('\t'.join(line_data) + '\n')
                                            
                    # add the information from the line to the summarized data dict for chromosome -- stats purposes
                    chrom = 'Scaf' if 'scaffold' in chrom.lower() else chrom
                    summarized_data[chrom]['sites'] += 1
                    summarized_data[chrom]['snps'] += int(snp_site['dp'][0])
                    summarized_data[chrom]['dp'] += snp_site['dp'][1]
                    summarized_data[chrom]['hets'] += snp_site['het']
                    summarized_data[chrom]['maf'][0] += snp_site['ref']
                    summarized_data[chrom]['maf'][1] += snp_site['alt']
                        
    summarized_data['totals'] = [total_snps, post_filter, population]

    return summarized_data, why
                        
                
def stat_filter(vcf, out_prefix, filter_params, post_filter_summary, stats_outfile, why_dropped):
    # make a stats sheet output for the filtering

    # write out the header info and filter parameters
    stats_outfile.write('AAFC VCF FILTER\nAuthor: Brian James (brian.james4@canada.ca)\nCreated: SEPT 2020\n\n')
    stats_outfile.write('Parameters as interpreted:\n\t')
    
    # filter_params are a dictionary, determine if False before writing line
    stats_outfile.write(f'--vcf {vcf.name}\n\t')
    stats_outfile.write(f'--out {out_prefix.name}\n\t')
    
    if filter_params['QUAL']:
        stats_outfile.write(f"--qual {filter_params['QUAL']}\n\t")
        
    if filter_params['MISS']:
        stats_outfile.write(f"--miss {filter_params['MISS']}\n\t")
        
    if filter_params['MAF']:
        stats_outfile.write(f"--maf {filter_params['MAF']}\n\t")
        
    if filter_params['MIN_DP']:
        stats_outfile.write(f"--mindp {filter_params['MIN_DP']}\n\t")
        
    if filter_params['DP_HOMO']:
        stats_outfile.write(f"--dphom {filter_params['DP_HOMO']}\n\t")
        
    if filter_params['DP_HET']:
        stats_outfile.write(f"--dphet {filter_params['DP_HET']}\n\t")
        
    if filter_params['CSV']:
        stats_outfile.write(f"--csv\n\t")
        
    if filter_params['INDEL']:
        stats_outfile.write(f"--removeindels\n\t")
    
    # write out the stats for before filter to post filter 
    stats_outfile.write(f"\n\nVCF FILTER STATS:\n")
    stats_outfile.write(f"After filtering, kept {post_filter_summary['totals'][1]} out of a possible {post_filter_summary['totals'][0]} sites\n")
    stats_outfile.write(f"Total Samples in VCF: {post_filter_summary['totals'][2]}\n\n")

    # include why lines are being dropped
    stats_outfile.write(f"SNP sites dropped due to the following parameters:\n")
    for reason, sites in why_dropped.items():
        stats_outfile.write(f'\t{reason}:\t{sites}\n')

    
    # write out the stats per chromosome
    stats_outfile.write('\nGenome Distribution of Remaining SNP Sites and Stats:\n\n')
    stats_outfile.write('CHROM\tSITES\t#SNPs\tAVG.DP\t%HET\tAVG.MAF\n')
    for key, val in post_filter_summary.items():
        try:
            avg_dp = float(val['dp'] / val['snps']) 
            phet =  float(val['hets'] / val['snps'])*100
            maf = float(1 - (val['maf'][0]/(val['maf'][0] + val['maf'][1])))
            
            stats_outfile.write(f"{key}\t{val['sites']}\t{val['snps']}\t{avg_dp:.2f}\t{phet:.2f}\t{maf:.3f}\n")
        except:
            continue

def main():
    args = get_arguments()
    vcf_file = Path(args.vcf) # initialize the file in the Path class
    outfile = args.outfile
    csv_out = args.csvout
    
    # create the outfile stuff required for the program, return a Path class for outfile and stats file
    vcf_out, stats_out = generate_outfile(vcf_file, outfile, csv_out)
    
    # add the filter parameters to a dictionary to pass into VCF parsing function
    filter_params = {'QUAL': int(args.min_quality), 'MISS': float(args.max_missing), 'MAF': float(args.minor_allele_freq), 
                     'INDEL': args.removeindels, 'DP_HOMO': int(args.min_depth_homozygous), 
                     'DP_HET': int(args.min_depth_heterozygous), 'MIN_DP': int(args.minimimum_snp_depth), 'CSV': csv_out}
    
    if csv_out:
        filter_params['INDEL'] = True
    
    # determine if VCF file is compressed or not, then filter the file using "filter_vcf"
    if '.gz' in vcf_file.name:
        with gzip.open(vcf_file, mode='rb', compresslevel=9, encoding=None, errors=None, newline=None) as vcf:
            if csv_out:
                with vcf_out.open('w', newline='') as outfile:
                    vcf_writer = csv.writer(outfile, delimiter=',')
                    stats_data, why_dropped = filter_vcf(vcf, filter_params, vcf_writer) 
            else:
                with vcf_out.open('w', newline=None) as outfile:
                    stats_data, why_dropped = filter_vcf(vcf, filter_params, outfile)

    else:
        with vcf_file.open('r') as vcf:
            if csv_out:
                with vcf_out.open('w', newline='') as outfile:
                    vcf_writer = csv.writer(outfile, delimiter=',')
                    stats_data = filter_vcf(vcf, filter_params, vcf_writer) 
            else:
                with vcf_out.open('w', newline=None) as outfile:
                    stats_data = filter_vcf(vcf, filter_params, outfile)
    
    # write the stats of the filtered VCF 
    with stats_out.open('w', newline=None) as outstats:
        stat_filter(vcf_file, outfile, filter_params, stats_data, outstats, why_dropped)


if __name__ == '__main__':
    # Steps for program:
    # 1. get the arguments for filtering
    # 2. make a new directory based on the output file name
    # 3. run the filter and output into new directory
    # 4. run the stats and output the stats file into new directory  
    
    main()
