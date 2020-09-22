from pathlib import Path
import argparse
import gzip
import numpy as np
from collections import defaultdict
import re
from math import floor
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib import ticker


def get_arguments():
    parser = argparse.ArgumentParser(description='Get a summary of vcf file depth')
    parser.add_argument('-v', '--vcf', required=True, help="Path to the VCF file.  VCF file can be compressed (gzip) or uncompressed")
    parser.add_argument('-w', '--window', required=True, help='Size of window for coverage analysis')
    args = parser.parse_args()
    return args


def windowed_variance(window, vcf_file):
    '''
    Go through the vcf file and return a sample depth per chromosome dictionary and total depth average
    :param vcf_file: Path Object for VCF file
    :param window:  Size of window to analyze the VCF file
    :return: Windowed Variance graph
    '''
    bin_size = int(window)
    zipped_file = True if '.gz' in vcf_file.name else False

    if zipped_file:
        with gzip.open(vcf_file, mode='rb', compresslevel=9, encoding=None, errors=None, newline=None) as vcf:
            
            window_data = {}
            for line in vcf:
                if line.decode('utf-8').startswith('##contig'):
                    # add the contig information to the dictionary (chromosome names and lengths)
                    if 'scaffold' not in line.decode('utf-8').lower():
                        chrom_name = re.findall(r'ID=(.*),', line.decode('utf-8'))[0]
                        chr_len = int(re.findall(r'length=(.*)>', line.decode('utf-8'))[0])
                        
                        # weird contigs in sunflower genome, maybe its a common thing to filter out?
                        if chrom_name == 'CP':
                            continue
                        elif chrom_name == 'MT':
                            continue
                        
                        else:
                            window_data[chrom_name] = {}
                            for p in range(bin_size, chr_len, bin_size):  
                                # add the bins from the window
                                window_data[chrom_name][p] = 0

                    else: # skip scaffolds
                        continue
                        
                elif line.decode('utf-8').startswith("#"):  
                    # header line, no info needed
                    continue           
                
                else:
                    snp = line.strip().decode('utf-8').split()

                    chrom = snp[0]
                    if 'scaffold' not in chrom.lower():
                    
                        # find the key in the dictionary to add this lines numbers to
                        pos = snp[1]
                        
                        # the nearest division of the window would be pos/window
                        
                        pos_key = int(floor(float(pos)/float(window))) * int(window)
                        if pos_key == 0:
                            pos_key = int(window)

                        # count the SNPs in the sample columns if called (ie not ./.)                       
                        for site in snp[9:]:  # just the sample information in a loop
                            snp_call = site.split(':')[0]
                            if snp_call != './.':
                                window_data[chrom][pos_key] += 1
                    else:
                        continue
    else:
        with vcf_file.open(mode='rb', newline=None) as vcf:
            
            window_data = {}
            for line in vcf:
                if line.decode('utf-8').startswith('##contig'):
                    # add the contig information to the dictionary (chromosome names and lengths)
                    if 'scaffold' not in line.decode('utf-8').lower():
                        chrom_name = re.findall(r'ID=(.*),', line.decode('utf-8'))[0]
                        chr_len = int(re.findall(r'length=(.*)>', line.decode('utf-8'))[0])
                        
                        # weird contigs in sunflower genome, maybe its a common thing to filter out?
                        if chrom_name == 'CP':
                            continue
                        elif chrom_name == 'MT':
                            continue
                        
                        else:
                            window_data[chrom_name] = {}
                            for p in range(bin_size, chr_len, bin_size):  
                                # add the bins from the window
                                window_data[chrom_name][p] = 0

                    else: # skip scaffolds
                        continue
                        
                elif line.decode('utf-8').startswith("#"):  
                    # header line, no info needed
                    continue           
                
                else:
                    snp = line.strip().decode('utf-8').split()

                    chrom = snp[0]
                    if 'scaffold' not in chrom.lower():
                    
                        # find the key in the dictionary to add this lines numbers to
                        pos = snp[1]
                        
                        # the nearest division of the window would be pos/window
                        
                        pos_key = int(floor(float(pos)/float(window))) * int(window)
                        if pos_key == 0:
                            pos_key = int(window)

                        # count the SNPs in the sample columns if called (ie not ./.)                       
                        for site in snp[9:]:  # just the sample information in a loop
                            snp_call = site.split(':')[0]
                            if snp_call != './.':
                                window_data[chrom][pos_key] += 1
                    else:
                        continue

    return window_data


def plot_windows(df, window):
 
    # dont plot the 0s, change them in the dataframe to nan
    df = df.replace(0, np.nan)
    
    chroms = list(df.columns)
    num_plots = len(chroms)
    cm = plt.get_cmap('viridis')  # fun color scheme
    color_list = [cm(1.*i/num_plots) for i in range(num_plots)] # generate the color list evenly across the color map viridis

    # get the X axis data split from scientific notation (just leave the value, remove the power) so we can add power to axis label instead
    # this just cleans the look of the data to not have notation in each subplot
    sci_not = '{:.2E}'.format(max(df.index)).split('E')[1] # just grab the power value ie) 2.0E+12 = +12.  This is used for X Axis Label

    x_axis_normalized = ['{:.1f}'.format(x/float(f'1E{sci_not}')) for x in df.index] 
    
    # get the number of subplots even to number of chromosomes
    fig, axes = plt.subplots(nrows=1, ncols=num_plots, sharey=True)
    fig.set_size_inches(12, 6, forward=True)
    fig.subplots_adjust(hspace=0, wspace=0)

    # set font for the title
    fig.suptitle(f'Genome Distribution of SNPs by Chromosome ({window}bp window)', fontsize='large', fontweight='bold')
    fig.text(0.5, 0.04, f'Position on Chromosome (1E{sci_not})', ha='center', fontsize='medium', fontweight='bold')
    fig.text(0.04, 0.5, "Number SNP's Per Window", va='center', rotation='vertical', fontsize='medium', fontweight='bold')

    for ax, y_data, chrom_color in zip(axes.flatten(), chroms, color_list):
        # pass in a series of zipped info to plot.  chroms is the name of the column, so to get Y data, call df[chrom] for each chrom in VCF
        ax.scatter(x_axis_normalized, df[y_data], color=chrom_color)
        ax.xaxis.set_major_locator(plt.MaxNLocator(2))
        ax.tick_params(axis='x', which='major', labelsize=6)
        ax.set_facecolor('aliceblue')
        ax.set_title(y_data.upper(), size=8, fontweight='bold')
    
    return plt


def main():
    args = get_arguments()
    vcf_file = args.vcf
    window = args.window

    vcf = Path(vcf_file)  # get the path to the VCF file into Path Class

    # generate the outfile save name (graph)
    if vcf.suffix == '.gz':
        vcf_with_suffix = vcf.stem
        just_vcf_name = Path(vcf_with_suffix).stem
    else:
        just_vcf_name = vcf.stem
    
    outdir = Path(vcf_file).parents[0]
    out = outdir / f'{just_vcf_name}.coverage.png'

    chrom_counts = windowed_variance(window, vcf)  # generates the windowed data (# SNPs per window)
    
    df = pd.DataFrame(chrom_counts) # convert to dataframe to visualize data
    coverage = plot_windows(df, window)

    # save the graph of coverage to the outfile directory.
    coverage.savefig(out)

    
if __name__ == '__main__':
    main()
