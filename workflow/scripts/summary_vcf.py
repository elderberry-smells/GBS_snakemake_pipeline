from pathlib import Path
import argparse


def get_arguments():
    parser = argparse.ArgumentParser(description='Get a summary of all VCFtools filter log files')
    parser.add_argument('-v', '--vcf_directory', help='path to the VCF filter directory where log files live')
    parser.add_argument('-o', '--out', help='Output of the summary')
    parser.add_argument('-s', '--split', help="True - VCF split into chromosomes, False - only one VCF")
    prog_args = parser.parse_args()
    return prog_args



def parse_vcflog(logfiles, summary):
    """
    Open the log file and get the information for a summary table
    :param logfile: the log file from vcftools
    :return: info from logfile
    """
    # format [plant, chrom, individuals start, individuals end, sites start, sites end]
    vcf_info = []
    sites = {"start": 0, "end": 0}
    for logfile in logfiles:
        filter_summary = []
        with open(logfile, 'r') as lfile:
            for line in lfile:
                if line.startswith("\t--vcf"):
                    data = line.split()
                    plant = data[1].split(".")[0]
                    chrom = data[1].split(".")[1]
                    filter_summary.append(plant)
                    filter_summary.append(chrom)
                elif line.startswith("After filtering"):
                    if "Individuals" in line:
                        data = line.split()
                        for i in data:
                            try:
                                indivs = int(i)
                                filter_summary.append(str(indivs))
                            except:
                                continue
                    else:
                        data = line.split()
                        nums = []
                        for i in data:
                            try:
                                nums.append(int(i))
                            except:
                                continue
                        filter_summary.append(str(nums[1]))
                        sites["start"] += nums[1]
                        filter_summary.append(str(nums[0]))
                        sites["end"] += nums[0]
        vcf_info.append(filter_summary)
    sites["retained"] = "{:.2f}".format((sites["end"]/sites["start"]) * 100)

    with open(summary, 'w') as vcf_sum:
        vcf_sum.write("VCF\tChr\tLines_start\tLines_end\tSites_start\tSites_end\n")
        for log_info in vcf_info:
            vcf_sum.writelines("\t".join(log_info))
            vcf_sum.write("\n")
        vcf_sum.write(f"\nTotal Sites Start:\t{sites['start']}\n"
                      f"Total Sites End:\t{sites['end']}\n"
                      f"Percent Retained:\t%{sites['retained']}")

def summarize_vcffilter(vcf_file, outfile):
    """
    Read the VCF file that was filtered and produce stats on the number of sites
    :param vcf_file: filtered vcf file
    :return: output file containing stats
    """
    vcf_name = Path(vcf_file[0])
    vcf_prefix = vcf_name.name.split('.')[0]


    site_count = {"total": 0}
    with open(vcf_file[0], 'r') as vcf:
        for line in vcf:
            if line.startswith('##') or line.startswith('#'):
                continue
            else:
                chr = line.split()[0]
                if chr in site_count:
                    site_count[chr] += 1
                    site_count["total"] += 1
                else:
                    site_count[chr] = 1
                    site_count["total"] += 1

    with open(outfile, "w") as vcf_sum:
        vcf_sum.write("VCF\tChr\tNumber_Sites\n")

        for key, val in site_count.items():
            if key != 'total':
                vcf_sum.write(f"{vcf_prefix}\t{key}\t{str(val)}\n")

        vcf_sum.write(f"\nTotal Sites:\t{site_count['total']}")


if __name__ == '__main__':
    # load arguments
    args = get_arguments()
    vcf_dir = args.vcf_directory
    outfile = args.out
    split_vcf = args.split

    if split_vcf.lower() == 'true':
        # get a list of all filter log files in directory
        log_files = list(Path(vcf_dir).glob('*.log'))
        parse_vcflog(log_files, outfile)

    elif split_vcf.lower() == 'false':
        vcf_file = list(Path(vcf_dir).glob('*.recode.vcf'))
        summarize_vcffilter(vcf_file, outfile)
    else:
        print("argument for split must be True or False")