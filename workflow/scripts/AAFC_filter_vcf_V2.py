#!/usr/bin/env python3
from pathlib import Path
from stats_output import summary_html
import csv
import re
import argparse
import gzip


def get_arguments():

    parser = argparse.ArgumentParser(description='Get a summary of vcf file depth')
    parser.add_argument('-vcf', '--vcf', required=True, help="VCF file you want to filter.  Can be compressed with gzip or uncompressed")
    parser.add_argument('-out', '--out_prefix', required=True, help="Prefix for filtered VCF file.  Will create a new directory with prefix name for output")
    parser.add_argument('-csv', '--csvout', action='store_true', help='After filtering, remove headers and return a human readable CSV version of the VCF')
    parser.add_argument('-qual', '--min_quality', default=False, help='Set the QUAL minimum for the SNP site, remove SNP if score is below [qual] argument')
    parser.add_argument('-miss', '--max_missing', default=False, help="Remove SNP site if site doesn't have AT LEAST [max_missing] calls across pop.")
    parser.add_argument('-maf', '--minor_allele_freq', default=False, help="Remove SNP site if minor allele frequency is lower than [minor_allele_freq]")
    parser.add_argument('-removeindels', '--removeindels', action='store_true', help="Remove SNP site if REF or ALT is an indel")
    parser.add_argument('-dphom', '--min_depth_homozygous', default=False, help="Do not count homozygous SNP for sample if depth doesn't match or exceed [min_depth_homozygous]")
    parser.add_argument('-dphet', '--min_depth_heterozygous', default=False, help="Do not count heterozygous SNP for sample if depth doesn't match or exceed [min_depth_homozygous]")
    parser.add_argument('-mindp', '--minimimum_snp_depth', default=False, help="Remove SNP site when snp depth is BELOW [minimimum_snp_depth], over_rides homo/het depth if higher")
    parser.add_argument('-biallelic', '--biallelic', action='store_true', help="Remove SNP site if the site not bi-allelic ie. max alleles 2 and min alleles 2 for REF/ALT")
    parser.add_argument('-summary', '--summary', action='store_true', help="Create an HTML summary output of the VCF file")
    args = parser.parse_args()
    return args


def generate_outfile(infile_path, out_prefix, csv_output):
    # take the arguments for infile/outfile and generate the directory and paths for the filter out and stats out
    # figure out here if the file or directory already exists, then inside the directory add the (1) or (2) or whatever copy this would be.
        
    # take the out argument and create a directory with that prefix in same directory as VCF file
    
    # get the path only from the VCF input Path class
    path_to_vcf = infile_path.parents[0]
    outdir = Path(path_to_vcf / f'{out_prefix}')
    outdir.mkdir(parents=True, exist_ok=True)
    
    # stats directory in the out directory
    stats_dir = Path(outdir / 'stats')
    stats_dir.mkdir(parents=True, exist_ok=True)
    
    # get just the name of the file, no suffix or filtered etc.
    vcf = infile_path.name.split('.')[0]  # use split on name in case of .gz or not
    
    if csv_output:
        vcf_out = Path(outdir / f'{vcf}.{out_prefix}.csv')
    else:
        vcf_out = Path(outdir / f'{vcf}.{out_prefix}.vcf')
    
    filter_stats = Path(outdir / f'{vcf}.{out_prefix}.log')
    
    stats_outpaths = {
        'chromsome_stats': Path(stats_dir / f'{vcf}.{out_prefix}.chr'),
        'snp_counts': Path(stats_dir / f'{vcf}.{out_prefix}.snps'),
        'DP':Path(stats_dir / f'{vcf}.{out_prefix}.depth'),
        'heterozygosity': Path(stats_dir / f'{vcf}.{out_prefix}.hets'),
        'allele_freq': Path(stats_dir / f'{vcf}.{out_prefix}.allele'),
        'summary': Path(outdir / f'{vcf}.{out_prefix}_summary.html')
    }    
    
    return vcf_out, filter_stats, stats_outpaths


def filter_vcf(infile, filter_params, outfiles):
  
    summarized_data = {}
    sample_stats = {}
    sample_list = []
    chromosomes = []
    population = 0
    total_snps = 0
    post_filter = 0
    
    why = {'QUAL': 0, 'INDEL': 0, 'BI-AL': 0, 'DP': 0, 'MISS': 0, 'MAF': 0}
    
    for line in infile:
        
        try:
            data = line.decode("utf-8").strip()  # decode because VCF was compressed and is in binary mode
        except:
            data = line.strip()  # just strip the line, no compresion in file
        
        if data.startswith('##'):
            if 'contig' in data.lower():  # isolates the header lines detailing the chromosomes and scaffold information
                
                if 'scaffold' not in data.lower():
                    skip_chroms = ['CP', 'MT']
                    chrom_name = re.findall(r'ID=(.*),', data)[0]
                    if chrom_name not in skip_chroms:
                        
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
                        if 'Scaf' in chromosomes:
                            continue
                        else:
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
            # make a sample dictionary to include a sample stats output, including population size
            for samp in data.split()[9:]:       
                sample_list.append(samp)
                population += 1  # add the sample # to the population list
                
            # add the samples to the sample_stats dictionary, and make a chromosome summary as well for final stats outputs
            counter = 9
            for samp in sample_list:
                sample_stats[counter] = {'sample':samp, 'stats':{}}
                for chrom in chromosomes:
                    # make a chromosome summary dictionary
                    summarized_data[chrom] = {'sites': 0, 'snps':0, 'dp': 0, 'hets': 0, 'maf': 0}
                    
                    # make a sample summary dictionary
                    if chrom != 'Scaf':
                        sample_stats[counter]['stats'][chrom] = {'snps':0, 'dp': 0, 'hets': 0, 'ref': 0, 'alt': 0, 'minor':0}
                        
                counter += 1
                              
            # write the line to the outfile
            if filter_params['CSV']:
                skip_col = ['INFO', 'FORMAT']
                try:  # make sure it is decoded if the file was gzipped
                    header = [i.decode("utf-8") for i in line.split() if i.decode("utf-8") not in skip_col]
                except:
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
                line_data = [snp[0], snp[1], snp[2], snp[3], snp[4], snp[5], snp[6]] 
            else:
                line_data = [snp[0], snp[1], snp[2], snp[3], snp[4], snp[5], snp[6], snp[7], snp[8]]  
            
            # vcf format for each line/site in file for easier reference
            chrom = snp[0]
            ref = snp[3]
            alt = snp[4]
            qual = float(snp[5])
            info = snp[7]
            biallelic = True if len(ref) == 1 and len(alt) == 1 else False
            
            # screen site based on INFO before moving into specific format info for each sample (if info set)
            if filter_params['QUAL'] and qual < int(filter_params['QUAL']):
                why['QUAL'] += 1
                continue
            
            elif filter_params['INDEL'] and 'INDEL' in info:
                why['INDEL'] += 1
                continue
            
            elif filter_params['BI'] and not biallelic:
                why['BI-AL'] += 1
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
                            line_data.append('N/A')
                        else:
                            line_data.append(site)

                    else:
                        if gt == './.':  # missing data point
                            snp_site['miss'] += 1
                            # add the info to the line to write(even though it may not be written in the end after filtering for MISS)
                            if filter_params['CSV']:
                                line_data.append('N/A')
                            else:
                                line_data.append(site)
                            
                        elif gt == '0/0':  # reference base called
                            if filter_params['DP_HOMO'] and int(dp) < int(filter_params['DP_HOMO']):
                                dp_drops['DP_HOMO'] += 1
                                site = site.replace('0/0:', './.:')  # didn't pass the test, make it no data and append it to snp_site for tally
                                snp_site['miss'] += 1
                                if filter_params['CSV']:
                                    line_data.append('N/A')
                                else:
                                    line_data.append(site)
                                snp_site['miss'] += 1
                                
                            else:  # passed the homozygous filter for ref
                                snp_site['ref'] += 2
                                snp_site['dp'][0] += 1  # add snp # tally
                                snp_site['dp'][1] += int(dp)  # add depth to make avg. depth at end
                                if filter_params['CSV']:
                                    line_data.append(f'{ref}:{ref}')
                                else:
                                    line_data.append(site)                        
                            
                        elif gt == '1/1':  # reference base called
                            if filter_params['DP_HOMO'] and int(dp) < int(filter_params['DP_HOMO']):
                                dp_drops['DP_HOMO'] += 1
                                site = site.replace('1/1:', './.:')  # didn't pass the test, make it no data and append it to snp_site for tally
                                if filter_params['CSV']:
                                    line_data.append('N/A')
                                else:
                                    line_data.append(site)
                                snp_site['miss'] += 1
                                               
                            else:  # passed the homozygous filter for ref
                                snp_site['alt'] += 2
                                snp_site['dp'][0] += 1  # add snp # tally
                                snp_site['dp'][1] += int(dp)  # add depth to make avg. depth at end
                                if filter_params['CSV']:
                                    line_data.append(f'{alt}:{alt}')
                                else:
                                    line_data.append(site)
                                                    
                        else:  # its a heterozygous SNP
                            if filter_params['DP_HET'] and int(dp) < int(filter_params['DP_HET']):
                                dp_drops['DP_HET'] += 1
                                site = site.replace('0/1:', './.:')  # didn't pass the test, make it no data and append it to snp_site for tally
                                
                                if filter_params['CSV']:
                                    line_data.append('N/A')
                                else:
                                    line_data.append(site)
                                    snp_site['miss'] += 1
                                                           
                            else:  # passed the heterozygous filter for ref
                                snp_site['het'] += 1
                                snp_site['dp'][0] += 1  # add snp # tally
                                snp_site['dp'][1] += int(dp)  # add depth to make avg. depth at end
                                
                                if filter_params['CSV']:
                                    line_data.append(f'{ref}:{alt}')
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
                    summarized_data[chrom]['maf'] += maf
                    
                    ####  SAMPLE INFORMATION STATS ####
                    # samples start on the 9th entry in each line of VCF, so calling the individual samples is done by looping from the 9th spot in line to end of list
                    for entry in range(9, len(sample_list)+9):
                        if chrom != 'Scaf':  # don't put the scaffold data into the sample info
                            sample_data = snp[entry]
                            gt = sample_data.split(':')[0]
                            dp = sample_data.split(':')[2]
                            
                            if gt != './.':  # collect data for sample is SNP was called
                                sample_stats[entry]['stats'][chrom]['snps'] += 1
                                sample_stats[entry]['stats'][chrom]['dp'] += int(dp)
                                if gt == '0/0':
                                    sample_stats[entry]['stats'][chrom]['ref'] += 1
                                    if snp_site['ref'] < snp_site['alt']:  
                                        # was the minor allele
                                        sample_stats[entry]['stats'][chrom]['minor'] += 1
                                elif gt == '1/1':
                                    sample_stats[entry]['stats'][chrom]['alt'] += 1
                                    if snp_site['ref'] > snp_site['alt']:  
                                        # was the minor allele
                                        sample_stats[entry]['stats'][chrom]['minor'] += 1
                                else:
                                    sample_stats[entry]['stats'][chrom]['hets'] += 1
                            else:
                                continue
                        else:
                            continue

    summarized_data['totals'] = [total_snps, post_filter, population]
    
    return summarized_data, why, sample_stats, chromosomes
                        
                
def filter_log(vcf, out_prefix, filter_params, post_filter_summary, stats_outfile, why_dropped):
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
        
    if filter_params['BI']:
        stats_outfile.write(f"--biallelic\n\t")
    
    # write out the stats for before filter to post filter 
    stats_outfile.write(f"\n\nVCF FILTER STATS:\n")
    stats_outfile.write(f"After filtering, kept {post_filter_summary['totals'][1]} out of a possible {post_filter_summary['totals'][0]} sites\n")
    stats_outfile.write(f"Total Samples in VCF: {post_filter_summary['totals'][2]}\n\n")

    # include why lines are being dropped
    stats_outfile.write(f"SNP sites dropped due to the following parameters:\n")
    for reason, sites in why_dropped.items():
        stats_outfile.write(f'{reason}:\t{sites}\n')


def filter_stats(vcf, stats_summary, sample_summary, chromosomes, stat_writers):
    # take the statistics from the VCF filter dictionaries, and write out the the files to include with the filtered VCF
       
    ### CHROMSOME STATISTICS FILE
    with stat_writers['chromsome_stats'].open('w', newline=None) as chr_writer:
        chr_writer.write('CHROM\tSITES\tSNPS\tAVGDP\tPHET\tAVGMAF\n')
        for key, val in stats_summary.items():
            # iterate through the keys (chromosome) and return stats on the values (counts, averages)
            try:
                avg_dp = float(val['dp'] / val['snps']) 
                phet =  float(val['hets'] / val['snps'])*100
                maf = val['maf']/val['sites']

                chr_writer.write(f"{key}\t{val['sites']}\t{val['snps']}\t{avg_dp:.2f}\t{phet:.2f}\t{maf:.3f}\n")

            except:  # avoids any contigs that have no data, but are still present in the VCF file
                continue


    ### THE FOLLOWING ARE FOR SAMPLE SPECIFIC STATISTICS, excluding the scaffolds

    sample_chroms = [i for i in chromosomes if i != 'Scaf'] # get sample stats header lines started
    
    ### SAMPLE - SNPS PER CHROMOSOME
    header = ['Sample']
    for chr in sample_chroms:
        
        header.append(chr)
    header.append('TotalSNPs')

    with stat_writers['snp_counts'].open('w', newline=None) as snp_writer:
        snp_writer.write('\t'.join(header) + '\n')
        
        for key, val in sample_summary.items():
            sample_snp_line = [val['sample']]

            total_snps = 0

            for x in sample_chroms:
                total_snps += val['stats'][x]['snps']
                    
                sample_snp_line.append(str(val['stats'][x]['snps']))
            
            sample_snp_line.append(str(total_snps))
            snp_writer.write('\t'.join(sample_snp_line) + '\n')
            
    ### SAMPLE - DP PER CHROMOSOME
    
    header = ['Sample']
    for chr in sample_chroms:
        header.append(chr)
    header.append('TotalAvgDP')

    with stat_writers['DP'].open('w', newline=None) as dp_writer:
        dp_writer.write('\t'.join(header) + '\n')
        
        for key, val in sample_summary.items():
            dp_line = [val['sample']]

            total_dp = [0, 0]

            for x in sample_chroms:
                total_dp[0] += val['stats'][x]['snps']
                total_dp[1] += val['stats'][x]['dp']
                
                num_snps = str(val['stats'][x]['snps'])
                
                if num_snps != '0':
                    avg_depth = val['stats'][x]['dp'] / val['stats'][x]['snps']
                else:
                    # no SNPs found, so don't get a divide by zero error.  This type of thing should be exception handling, but thats for later
                    avg_depth = 0
                    
                dp_line.append(f'{avg_depth:.1f}')
                
            if total_dp[0] != 0:
                percent_dp = (total_dp[1] / total_dp[0])
            else:  # your filter reduced line to no samples, this avoids divide by 0 error
                percent_dp = 0
                
            dp_line.append(f'{percent_dp:.1f}')
            dp_writer.write('\t'.join(dp_line) + '\n')
    
    ### SAMPLE - HETEROZYGOSITY
    header = ['Sample']
    for chr in sample_chroms:
        header.append(f'{chr}')
    header.append('TotalPHet')
    
    with stat_writers['heterozygosity'].open('w', newline=None) as het_writer:
        het_writer.write('\t'.join(header) + '\n')
        
        for key, val in sample_summary.items():
            sample_het_line = [val['sample']]

            het_tally = [0, 0]  # tally of total snps, total het calls

            for x in sample_chroms:
                if val['stats'][x]['snps']!= 0:
                    percent_het = (val['stats'][x]['hets'] / val['stats'][x]['snps']) * 100
                else:
                    percent_het = 0
                
                sample_het_line.append(f'{percent_het:.1f}')
                # add info to tally for total % calculation
                het_tally[0] += val['stats'][x]['snps']
                het_tally[1] += val['stats'][x]['hets']
                
            if het_tally[0] != 0:
                total_het = (het_tally[1] / het_tally[0]) * 100
            else: # your filter reduced number to 0, this avoids divison by 0 error
                total_het = 0
                
            sample_het_line.append(f'{total_het:.1f}')
            
            het_writer.write('\t'.join(sample_het_line) + '\n')
            
    ### SAMPLE - Allele Freq
    header = ['Sample']
    for chr in sample_chroms:
        if 'Scaf' not in chr:
            header.append(f'{chr}')
    header.append('TotalMAF')
    
    with stat_writers['allele_freq'].open('w', newline=None) as allele_writer:
        allele_writer.write('\t'.join(header) + '\n')
        
        for key, val in sample_summary.items():
            allele_line = [val['sample']]

            minor_tally = [0, 0]  # tally of total snps, total minor allele calls

            for x in sample_chroms:
                allele_line.append(str(val['stats'][x]['minor']))

                # add info to tally for total % calculation
                minor_tally[0] += val['stats'][x]['snps']
                minor_tally[1] += val['stats'][x]['minor']
            
            if minor_tally[0] != 0:
                percent_minor = minor_tally[1] / minor_tally[0]
            else: # your filter reduced number to zero, this avoids division by zero error
                percent_minor = 0
                
            allele_line.append(f'{percent_minor:.2f}')
            
            allele_writer.write('\t'.join(allele_line) + '\n')
            
            
def main():
    args = get_arguments()
    vcf_file = Path(args.vcf) # initialize the file in the Path class
    outfile = args.out_prefix
    csv_out = args.csvout
    
    # create the outfile stuff required for the program, return a Path class for outfile and stats file
    vcf_out, filterlog_out, stats_out = generate_outfile(vcf_file, outfile, csv_out)
    
    # add the filter parameters to a dictionary to pass into VCF parsing function
    filter_params = {'QUAL': int(args.min_quality), 'MISS': float(args.max_missing), 'MAF': float(args.minor_allele_freq), 
                     'INDEL': args.removeindels, 'DP_HOMO': int(args.min_depth_homozygous), 'BI': args.biallelic,
                     'DP_HET': int(args.min_depth_heterozygous), 'MIN_DP': int(args.minimimum_snp_depth), 'CSV': csv_out,
                     'SUMMARY': args.summary}
    
    # if you want to keep indels and SNPs that aren't biallelic, comment out the next 3 lines with # 
    if csv_out:
        filter_params['INDEL'] = True
        filter_params['BI'] = True
    
    # determine if VCF file is compressed or not, then filter the file using "filter_vcf"
    if '.gz' in vcf_file.name:
        with gzip.open(vcf_file, mode='rb', compresslevel=9, encoding=None, errors=None, newline=None) as vcf:
            if csv_out:
                with vcf_out.open('w', newline='') as outfile:
                    vcf_writer = csv.writer(outfile, delimiter=',')
                    stats_data, why_dropped, sample_statistics, chromosomes = filter_vcf(vcf, filter_params, vcf_writer) 
            else:
                with vcf_out.open('w', newline=None) as outfile:
                    stats_data, why_dropped, sample_statistics, chromosomes = filter_vcf(vcf, filter_params, outfile)

    else:
        with vcf_file.open('r') as vcf:
            if csv_out:
                with vcf_out.open('w', newline='') as outfile:
                    vcf_writer = csv.writer(outfile, delimiter=',')
                    stats_data, why_dropped, sample_statistics, chromosomes = filter_vcf(vcf, filter_params, vcf_writer) 
            else:
                with vcf_out.open('w', newline=None) as outfile:
                    stats_data, why_dropped, sample_statistics, chromosomes = filter_vcf(vcf, filter_params, outfile)
    
    # write the log file for the filter
    with filterlog_out.open('w', newline=None) as outlog:
        filter_log(vcf_file, outfile, filter_params, stats_data, outlog, why_dropped)
        
    # write the statistics for the filtered VCF file from the numerous files you need to write
    filter_stats(vcf_file, stats_data, sample_statistics, chromosomes, stats_out)

    if filter_params['SUMMARY']:
        # things needed for summary:
        # 1. vcf file name 
        # 2. reasons why the lines were dropped (why dictionary)
        # 3. outfiles in stats - dictionary where you can grab the file names, including HTML output name (stats_out for that)
        summary_html(vcf_file, why_dropped, stats_data, stats_out)

if __name__ == '__main__':
    # Steps for program:
    # 1. get the arguments for filtering
    # 2. make a new directory based on the output file name
    # 3. run the filter and output into new directory
    # 4. run the stats and output the stats file into new directory  
    
    main()
