#!usr/bin/env python3
"""
This script is used to demultiplex paired end fastq format GBS files into the seperate GBS index
files based on the in-line barcodes annealed to the DNA fragments in Read1 files only

(The GBS data produced is accomplished by following Poland et. al 2012 GBS protocol, and sequenced on Illumina platform)

To run this script individually and not part of a pipeline use the following example command:

$ python3 scripts/PE_fastq_demultiplex.py -f path/to/read_1.fastq -s path/to/samplesheet.txt -b path/to/barcode_file.txt

Author:  Brian James
email:  brian.james4@canada.ca
"""

import re
import argparse
import os
from datetime import datetime
import concurrent.futures
import multiprocessing
from multiprocessing import Lock


def get_arguments():
    parser = argparse.ArgumentParser(description='Demultiplex a fastq file into individual sample files')
    parser.add_argument('-f', '--fastq', help='fastq file. Format: path/to/file.fastq.  Use R1 path')
    parser.add_argument('-b', '--barcode', help='Barcode file.  Format: path/to/barcode_file.txt')
    parser.add_argument('-s', '--sample', help='Sample config file.  Format: path/to/sample_file.txt')
    prog_args = parser.parse_args()
    return prog_args


def make_dicts(barcode_arg, sample_arg):
    """
    Taking in the barcode file and the sample file to create dictionaries for those values
    :param barcode_arg: barcode.txt file from arguments housing the index names and barcodes
    :param sample_arg: samplesheet.txt file from arguments housing sample name and index name
    :return: barcode_dict - dict format:  Index_name : [index_barcode, counter]
             sample_dict - dict format:  Index_name : sample_name
    """
    # get barcode and sample file into a dictionary format for passing into parse_barcode() function
    barcode_dict = {}
    with open(barcode_arg, 'r') as bar:
        content = bar.readlines()

        for line in content:
            data = line.split()
            barcode_dict[data[0]] = [data[1], 0]  # dict format:  Index_name : [index_barcode, counter]
            barcode_dict['unmatched'] = ['unmatched', 0]  # add unmatched to the dict (not in barcode sheet)

    # sample dictionary is formatted as 4 columns (tab delimited).  sample_number -- Index_name -- sample_name -- taxon
    sample_dict = {}

    with open(sample_arg, 'r') as samps:
        next(samps)  # skip the header line
        index_sample = samps.readlines()

        for line in index_sample:
            sam = line.split()
            if len(sam) == 0:  # avoid any blank lines in the sample sheet that might be hiding at end of file.
                continue
            else:
                sample_dict[sam[1]] = [sam[2], []]
        sample_dict['unmatched'] = ['unmatched', []]

    return barcode_dict, sample_dict


def process_fq(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}


def fq_write_lock(fq_writelock, read1, read2):
    """
    Lock the files you are writing to to avoid any Non-ThreadSafe writing to the file when you have a lot of processes
    accessing the same fiel to append to.
    :param lock: multiprocessing lock
    :param read1: sample dict containing the write data for parsed read 1 data
    :param read2: sample dict containing the write data for parsed read 2 data
    :return: no return, just write to the files.
    """
    fq_writelock.acquire()

    for samp, seq in read1.items():
        write_name = f'{path_name}/demultiplex/{seq[0]}.1.fastq'

        with open(write_name, 'a') as f1_out:
            f1_out.write(''.join(seq[1]))

    for samp2, seq2 in read2.items():
        write_name2 = f'{path_name}/demultiplex/{seq2[0]}.2.fastq'

        with open(write_name2, 'a') as f2_out:
            f2_out.write(''.join(seq2[1]))

    fq_writelock.release()


def parse_index(seq_record, bcode_dict):
    """
    :param seq_record: the record parsed from SeqIO.parse, houses info from fastq file
    :param bcode_dict:  dictionary of the information on the barcode/index name (index001: bcode, index2: bcode, ...)
    :return: returns the filename/index name from matching up index with barcode file
    """
    
    longest = max(len(bcode[0]) for bindex, bcode in bcode_dict.items())
    shortest = min(len(bcode[0]) for bindex, bcode in bcode_dict.items())

    sequence = str(seq_record['sequence'])
    bcode = re.search('TGCA', sequence)  # get everything up to (and including) the TGCA of the barcode

    if bcode:
        numbers = bcode.span()  # get the start and finish locations of TGCA

        # some barcodes have 2 TGCA in sequence, need a specific number of bases before TGCA to be proper depending on bcode_file.
        if longest == 14:   # using the newer set of barcodes.  Need a minimum of 5 bases before TGCA to properly match!
            pre_tgca = 5
        else:               # using the older set of barcodes.  Need a minimum of 4 bases before TGCA to properly match!
            pre_tgca = 4
            
        if numbers[0] < pre_tgca:  # not the proper amount of bases before TGCA found.  Check if another TGCA in sequence.
            bcode = re.search('TGCA', sequence[numbers[1]:])
            if bcode:
                new_span = bcode.span()
                numbers = (new_span[0] + numbers[1], new_span[1] + numbers[1])
            else:
                return 'unmatched', 0

        if numbers[1] in range(shortest, longest+1):  # make sure the result is the proper length of the barcodes provided

            index = sequence[0:numbers[1]]  # grab just the index + TGCA from sequence
            cut_index = numbers[0]
            bcode_len = len(index)

            if index.startswith('N'):
                # HiSeq will sometimes output a sequence that starts with N placeholder (instead of an ACGT)
                # we want to account for that one error and not immediately throw the data out.

                # get all partial matches for the barcode (excluding the first letter) ex) {'gbsx001':[TGCCCTATGCA, 0]}
                gbs_index = dict(filter(lambda item: index[1:] in item[1][0], bcode_dict.items()))

                # only keep the matches that are the same length as the index (if multiple hits, reduce to 1)
                for key, val in gbs_index.items():
                    if len(val[0]) == bcode_len:
                        return key, cut_index

                return 'unmatched', 0  # if still no match for length, return unmatched

            else:
                # get all partial matches for the barcode (including the first letter)
                gbs_index = dict(filter(lambda item: index in item[1][0], bcode_dict.items()))

                # only keep the matches that are the same length as the index (if multiple hits, reduce to 1)
                for key, val in gbs_index.items():
                    if len(val[0]) == bcode_len:
                        return key, cut_index

                return 'unmatched', 0  # if still no match for length, return unmatched

        else:  # contains TGCA, but can't determine based on size of index returned (not between 9 and 14 bp long)
            return 'unmatched', 0

    else:  # no TGCA was found, so no barcode can be matched.
        return 'unmatched', 0


def split_fq(fastq):

    fastq2 = fastq.replace('_R1', '_R2')  # get read 2

    # make the barcode and sample dictionaries
    barcodes, samples = make_dicts(barcode_file, sample_sheet)
    _, samples2 = make_dicts(barcode_file, sample_sheet)

    with open(fastq) as f1, open(fastq2) as f2:

        lines1 = []  # empty list to append one fastq entry to from read 1
        lines2 = []  # empty list to append one fastq entry to from read 2

        #  process read 1 entry and assign an index
        for r1, r2 in zip(f1, f2):
            lines1.append(r1.rstrip())
            lines2.append(r2.rstrip())

            if len(lines1) == 4:
                # get both records of read 1 and read 2:  returns dict ['name', 'sequence', 'optional', 'quality']
                record = process_fq(lines1)
                record2 = process_fq(lines2)

                # pass the records through the parse_fq function to get the index name for read 1
                index_name, slicer = parse_index(record, barcodes)
                barcodes[index_name][1] += 1

                # read 1 sliced data
                sequence = str(record['sequence'][slicer:])  # cut off the index, leaving just TCGA at start
                quality = str(record['quality'][slicer:])  # cut off the same amount of quality scores

                # read 2 sliced data
                sequence2 = str(record2['sequence'][slicer:])  # cut off the index, leaving just TCGA at start
                quality2 = str(record2['quality'][slicer:])  # cut off the same amount of quality scores

                # double check that the run_ids match between read1 and read 2
                run_id = record['name'].split(' ')[0]
                run_id2 = record2['name'].split(' ')[0]

                if run_id == run_id2:
                    # add the record information to the samples dictionary as a string
                    fastq1 = f"{record['name']}\n{sequence}\n+\n{quality}\n"
                    fastq2 = f"{record2['name']}\n{sequence2}\n+\n{quality2}\n"

                    # sometimes they match to index not actually in the sample sheet, so weed those ou
                    if index_name in samples:
                        # add them to the dictionaries value (list)
                        samples[index_name][1].append(fastq1)
                        samples2[index_name][1].append(fastq2)

                    else:  # assign it to unmatched instead
                        samples['unmatched'][1].append(fastq1)
                        samples2['unmatched'][1].append(fastq2)

                    lines1 = []
                    lines2 = []

                else:
                    print('your fastq file formats are not matching.  could be issues with fastq files')
                    continue

    # write the files for the sample dictionary entries for read 1
    fq_write_lock(lock, read1=samples, read2=samples2)

    return barcodes


def write_outstats(st_time, fq_path, bcode_file, bcode_dict):
    """
    :param st_time:  time that the program started
    :param fq_path:  The -f argument call, which includes the path in which to create and save the log
    :param bcode_file:  the -b argument call.  Add to the header info on what was run
    :param bcode_dict:  Barcode dictionary made from barcode file, includes barcode and counter as values for keys
    :return:  writes/appends the command, start and finish, and final output of demultiplex script to log file
    """

    if os.path.isdir(f'{os.path.split(fq_path)[0]}/log'):
        pass
    else:
        os.mkdir(f'{os.path.split(fq_path)[0]}/log')

    write_name = f'{os.path.split(fq_path)[0]}/log/demultiplex.output.log'

    command = f'Working Directory:  {os.getcwd()}\n'\
              f'Running Command:\n\t$ python3 scripts/fastq_demultiplex.py -f {fq_path} -b {bcode_file}\n'

    if st_time is not None:
        with open(write_name, 'w') as outfile:
            outfile.write(command)
            outfile.write(f'Started demultiplexing @: {st_time}\n')

    else:
        with open(write_name, 'a') as outfile:
            outfile.write(f'Finished demultiplexing @: {datetime.now()}\n\n')

            total_rcount = 0
            total_pcount = 1
            for bcode_key, bcode_val in bcode_dict.items():
                if bcode_key == 'unmatched':
                    total_rcount += bcode_val[1]
                else:
                    write_stats = f'{bcode_key}\t{bcode_val[0]}\t{bcode_val[1]}\n'
                    outfile.write(write_stats)
                    total_pcount += bcode_val[1]
                    total_rcount += bcode_val[1]

            percent_partitioned = (total_pcount / total_rcount) * 100
            reads = f'\nTotal reads in PE files: {total_rcount}\n'
            partitioned = f'Total partitioned reads: {total_pcount} \n'
            percent_p = f'Percentage partitioned : {percent_partitioned:.2f}'
            outfile.write(reads)
            outfile.write(partitioned)
            outfile.write(percent_p)


if __name__ == '__main__':

    # load the arguments from the command line as the functions going forward
    args = get_arguments()
    fastq_file1 = args.fastq                        # read 1.fastq file name
    barcode_file = args.barcode                     # barcode file (tab delimited, 2 headers)
    sample_sheet = args.sample                      # sample config file (tab delimited, 3 headers)
    start_time = datetime.now()                     # get the start time for the program starting
    available_cpus = multiprocessing.cpu_count()    # get the number of CPUs available for multiprocessing
    lock = Lock()

    # write header and info for stats output/log file
    write_outstats(st_time=start_time, fq_path=fastq_file1, bcode_file=barcode_file, bcode_dict=None)

    # check to see if demultiplex folder exists for this sample, if not, create it
    path_name, _ = os.path.split(fastq_file1)
    if os.path.isdir(f'{path_name}/demultiplex/'):
        pass
    else:
        os.mkdir(f'{path_name}/demultiplex/')

    demultiplexed_stats, _ = make_dicts(barcode_file, sample_sheet)

    # get a list of chunked data to demultiplex(just read 1 is fine)
    fq_chunks = os.listdir(f'{path_name}/chunks')
    fq_read1 = [f'{path_name}/chunks/{i}' for i in fq_chunks if '_R1' in i]

    ############### demultiplex read1 and read2 fastq data using multiprocessing tool concurrent.futures  ##############

    with concurrent.futures.ProcessPoolExecutor(max_workers=available_cpus) as executor:
        # pass in list of fastqs to process and be assigned to other CPUs
        future_to_fq = executor.map(split_fq, fq_read1)

        for result in future_to_fq:
            # update the demultiplexed stats
            for kb, vb in result.items():
                demultiplexed_stats[kb][1] += int(vb[1])

    ####################################  write out stats, close program ###############################################

    # write the stdout for the barcodes founds in the file.  total reads, etc.
    write_outstats(st_time=None, fq_path=fastq_file1, bcode_file=barcode_file, bcode_dict=demultiplexed_stats)

    # a try except clause to say that R1 is not in name could be helpful here?
