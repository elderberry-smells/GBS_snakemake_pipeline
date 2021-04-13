#!/usr/bin/env python3
from __future__ import print_function
from pathlib import Path
import pandas as pd
import numpy as np
import argparse
from jinja2 import Environment, FileSystemLoader
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff

def summary_html(vcf_file, why_dropped, stats_data, stats_out):
    '''
    Take in the info from the filter and produce a HTML representation of the data
    params vcf_file: the pathlib object detailing the info for the VCF file to be filtered
    params why_dropped: a dictionary that details the reasons why a SNP site did not pass the filter criteria
    params stats_data: a dictionary that houses the information like number of samples, pre filter and post filter sites 
    params stats_out: a dictionary that has all the file names in pathlib objects to grab summary data, and the HTML output path object
    return: nothing, just write the HTML file summary into same directory as VCF filter
    '''

    # filter summary
    samps = stats_data['totals'][2]
    pre_filt = stats_data['totals'][0]
    post_filt = stats_data['totals'][1]
    percent_retained = (post_filt/pre_filt) * 100
    pretain = f'{percent_retained:.2f}'

    #################### stacked plots for genome summary  ##################################################
    chr_df = pd.read_csv(stats_out['chromsome_stats'], delimiter='\t')

    genome_fig = make_subplots(rows=5, cols=1, shared_xaxes=True,
                    subplot_titles=("SNP Sites Per Chrom",
                                    "SNP Density Per Chrom", 
                                    "Average DP per Chrom", 
                                    "Average HET per chrom", 
                                    "Average MAF per chrom"))

    genome_fig.append_trace(go.Bar(
        x=chr_df['CHROM'],
        y=chr_df['SITES'],
        name='Sites'), row=1, col=1)

    genome_fig.append_trace(go.Bar(
        x=chr_df['CHROM'],
        y=chr_df['SNPS'],
        name='SNP Density'), row=2, col=1)

    genome_fig.append_trace(go.Scatter(
        x=chr_df['CHROM'],
        y=chr_df['AVGDP'],
        name='DP'), row=3, col=1)

    genome_fig.append_trace(go.Scatter(
        x=chr_df['CHROM'],
        y=chr_df['PHET'],
        name='Heterozygosity(%)'), row=4, col=1)

    genome_fig.append_trace(go.Scatter(
        x=chr_df['CHROM'],
        y=chr_df['AVGMAF'],
        name='MAF'), row=5, col=1)

    genome_fig.update_layout(height=900, width=900, title_text="Genome Summary of VCF")


    ############################## heterozygosity stuff  #####################################################
    het_df = pd.read_csv(stats_out['heterozygosity'], delimiter='\t')

    # get top 5 and bottom 5 het % in stats
    top5 = het_df.sort_values(by=['TotalPHet'], ascending=False).head(5)
    top5 = top5[['Sample', 'TotalPHet']].copy()
    bot5 = het_df.sort_values(by=['TotalPHet'], ascending=True).head(5)
    bot5 = bot5[['Sample', 'TotalPHet']].copy()
    mean_het = het_df['TotalPHet'].mean()
    het_cols = list(het_df.columns)
    het_cols = [i for i in het_cols if i not in ['Sample', 'Scaf']]

    # make a dist plot of heterozygosity that shows the distribution per chromsome and total %, where you can click the group and see that dist in the front.
    # done in bin sizes of 2%.
    hist_data = []
    group_labels = []
    for col in het_cols:
        het_data = het_df[col]
        hist_data.append(het_data)
        group_labels.append(col)
            
    het_fig = ff.create_distplot(hist_data, group_labels, curve_type="normal", bin_size=2, show_rug=False)
    het_fig.update_layout(title_text="Distribution of Samples based on % Heterozygosity") # Set title


    ################################  DEPTH STUFF ##################################################
    dp_df = pd.read_csv(stats_out['DP'], delimiter='\t')

    # get top 5 and bottom 5 DP samples in stats
    dp_top5 = dp_df.sort_values(by=['TotalAvgDP'], ascending=False).head(5)
    dp_top5 = dp_top5[['Sample', 'TotalAvgDP']].copy()
    dp_bot5 = dp_df.sort_values(by=['TotalAvgDP'], ascending=True).head(5)
    dp_bot5 = dp_bot5[['Sample', 'TotalAvgDP']].copy()
    mean_dp = dp_df['TotalAvgDP'].mean()

    # make a scatter plot of the depth data (each sample represented) and be able to select data based on chromosome or overall
    dp_fig = go.Figure()
    dp_fig.add_trace(go.Bar(x=dp_df['Sample'], y=dp_df['TotalAvgDP']))
    dp_fig.update_layout(title='Average Depth per Sample (across all chromosomes)',
                    yaxis_zeroline=False, xaxis_zeroline=False)

    ################################  MAF  #########################################################
    maf_df = pd.read_csv(stats_out['allele_freq'], delimiter='\t')

    maf_df['MAF'] = maf_df['TotalMAF'] * 100

    # get top 5 and bottom 5 MAF in stats
    maf_top5 = maf_df.sort_values(by=['TotalMAF'], ascending=False).head(5)
    maf_top5 = maf_top5[['Sample', 'TotalMAF']].copy()
    maf_bot5 = maf_df.sort_values(by=['TotalMAF'], ascending=True).head(5)
    maf_bot5 = maf_bot5[['Sample', 'TotalMAF']].copy()
    mean_maf = maf_df['TotalMAF'].mean()

    maf_hist_data= [maf_df['MAF']]
    maf_group_labels = ['MAF distribution']
    maf_fig = ff.create_distplot(maf_hist_data, maf_group_labels, show_rug=False)
    maf_fig.update_layout(title_text="Distribution of Samples based on Minor Alelle Frequency (converted to %)", ) # Set title

    ###  initialize the HTML template, fill in the variables, and write out the HTML file 
    env = Environment(loader=FileSystemLoader('/home/bioinf/vcf/scripts'))
    template = env.get_template("base2.html") # when added to cluster, put in full path
    template_vars = {"vcf_file": vcf_file.name,
                     "num_samp": samps,
                    "chr_summ" : chr_df.to_html(index=False, col_space=20, border=1, justify='center', table_id='chr'),
                    "graph": genome_fig.to_html(full_html=False, include_plotlyjs='cdn'),
                    "why_dropped": why_dropped,
                    "presites": pre_filt,
                    "postsites": post_filt,
                    "pretain": pretain,
                    "mean_dp": f'{mean_dp:.2f}',
                    "dp_top5": dp_top5.to_html(index=False, col_space=20, border=1, justify='center', table_id='dptop5'),
                    "dp_bot5": dp_bot5.to_html(index=False, col_space=20, border=1, justify='center', table_id='dpbot5'),
                    "dp_dist": dp_fig.to_html(full_html=False, include_plotlyjs='cdn'),
                    "mean_het": f'{mean_het:.2f}',
                    "het_top5": top5.to_html(index=False, col_space=20, border=1, justify='center', table_id='hettop5'),
                    "het_bot5": bot5.to_html(index=False, col_space=20, border=1, justify='center', table_id='hetbot5'),
                    "het_dist": het_fig.to_html(full_html=False, include_plotlyjs='cdn'),
                    "mean_maf": f'{mean_maf:.2f}',
                    "maf_top5": maf_top5.to_html(index=False, col_space=20, border=1, justify='center', table_id='maftop5'),
                    "maf_bot5": maf_bot5.to_html(index=False, col_space=20, border=1, justify='center', table_id='mafbot5'),
                    "maf_dist": maf_fig.to_html(full_html=False, include_plotlyjs='cdn')
                    }

    # render the template into HTML format with all the information generated above, then write it out to .html file
    html_out = template.render(template_vars)

    with stats_out['summary'].open('w') as f:
        f.write(html_out)
        