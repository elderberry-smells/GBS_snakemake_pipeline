#!usr/bin/env python3
"""
This script is used to demultiplex paired end fastq format GBS files into the seperate GBS index
files based on the in-line barcodes annealed to the DNA fragments in Read1 files only

To run this script individually and not part of a pipeline use the following example command:

$ python3 PE_fastq_demultiplex_AAFC.py -f read_1.fastq -s samplesheet.csv -b barcode_file.txt

Author:  Brian James
email:  brian.james4@canada.ca
"""

import re
import argparse
import os
from datetime import datetime
from random import randint
import concurrent.futures
import multiprocessing
from multiprocessing import Lock
import csv
from pathlib import Path


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
             sample_dict - dict format:  {Index_name : [ reference_dir_name, sample_ID, [empty_list] ]}
    """
    # get barcode and sample file into a dictionary format for passing into parse_barcode() function
    barcode_dict = {}
    with open(barcode_arg, 'r') as bar:
        content = bar.readlines()

        for line in content:
            data = line.split()
            barcode_dict[data[0]] = [data[1], 0]  # dict format:  Index_name : [index_barcode, counter]
            barcode_dict['unmatched'] = ['unmatched', 0]  # add unmatched to the dict (not in barcode sheet)

    # sample sheet is formatted as 4 columns (tab or comma delimited).
    # [0]sample_number -- [1]Index_name -- [2]sample_name -- [3]reference_file
    sample_dict = {}

    s_path, s_name = os.path.split(sample_arg)
    samp_ext = os.path.splitext(s_name)[1]

    # open the samplesheet, determine the format to read the data, and add the information to the sample_dict
    # sample dict format per sample:  {Index_name : [ reference_dir_name, sample_ID, [empty_list] ]}

    with open(sample_arg, newline='') as samps:
        next(samps)  # skip the header line
        if samp_ext == '.csv':
            reader = csv.reader(samps, delimiter=',')
        else:
            reader = csv.reader(samps, delimiter='\t')

        for row in reader:
            # get the information from the reference column, if empty, return unreferenced as directory
            ref_path, _ = os.path.split(row[3].strip())
            if ref_path != '':
                ref_dir = ref_path.split('/')[-1]
            else:
                ref_dir = 'unreferenced'

            sample_dict[row[1]] = [ref_dir, row[2], []]

        sample_dict['unmatched'] = ['unmatched', 'unmatched', []]

    return barcode_dict, sample_dict


def process_fq(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}


def fq_write_lock(fq_writelock, read1, read2):
    """
    Lock the files you are writing to to avoid any Non-ThreadSafe writing to the file when you have a lot of processes
    accessing the same fiel to append to.
    :param fq_writelock: multiprocessing lock
    :param read1: sample dict containing the write data for parsed read 1 data
    :param read2: sample dict containing the write data for parsed read 2 data
    :return: no return, just write to the files.
    """
    fq_writelock.acquire()

    for samp, seq in read1.items():
        write_name = f"{path_name}/{seq[0]}/demultiplex/{seq[1]}.1.fastq"

        with open(write_name, 'a') as f1_out:
            f1_out.write(''.join(seq[2]))

    for samp2, seq2 in read2.items():
        write_name2 = f'{path_name}/{seq2[0]}/demultiplex/{seq2[1]}.2.fastq'

        with open(write_name2, 'a') as f2_out:
            f2_out.write(''.join(seq2[2]))

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

        # some barcodes have 2 TGCA in sequence, need a specific number of bases before TGCA to be proper
        if longest == 14:   # using the newer set of barcodes.  Need a minimum of 5 bases before TGCA to properly match!
            pre_tgca = 5
        else:               # using the older set of barcodes.  Need a minimum of 4 bases before TGCA to properly match!
            pre_tgca = 4
            
        if numbers[0] < pre_tgca:  # not the proper amount of bases before TGCA found. Check if another TGCA in sequence
            bcode = re.search('TGCA', sequence[numbers[1]:])
            if bcode:
                new_span = bcode.span()
                numbers = (new_span[0] + numbers[1], new_span[1] + numbers[1])
            else:
                return 'unmatched', 0

        if numbers[1] in range(shortest, longest+1):  # make sure the barcode result is the proper length

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


def parse_riptide(seq_record, bcode_dict):
    """
    :param seq_record: the record is the sequence entry (4 lines) from fastq file
    :param bcode_dict:  dictionary of the information on the barcode/index name (index001: bcode, index2: bcode, ...)
    :return: returns the filename/index name from matching up index with barcode file
    """

    sequence = str(seq_record['sequence'])
    bcode = sequence[:8]

    return_name = []

    if bcode.startswith('N'):
        gbs_index = dict(filter(lambda item: bcode[1:] in item[1][0], bcode_dict.items()))

        # only keep the matches that are the same length as the index (if multiple hits, reduce to 1)
        for key, val in gbs_index.items():
            return_name.append(key)
    else:
        gbs_index = dict(filter(lambda item: bcode in item[1][0], bcode_dict.items()))

        for key, val in gbs_index.items():
            return_name.append(key)

    if len(return_name) == 0:
        return 'unmatched'
    else:
        return return_name[0]


def split_fq(fastq):

    fastq2 = fastq.replace('_R1', '_R2')  # get read 2

    # make the barcode and sample dictionaries
    barcodes, samples = make_dicts(barcode_file, sample_sheet)
    _, samples2 = make_dicts(barcode_file, sample_sheet)
    rand_x = int(randint(1000000, 9000000))  # random x_coord to start on, so each multiprocess doesn't make same numbers
    rand_y = int(randint(1000000, 9000000))  # random y_coord to start on, so each multiprocess doesn't make same numbers

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
                if 'riptide' in barcode_file:
                    index_name = parse_riptide(record, barcodes)
                    barcodes[index_name][1] += 1

                    # read 1 sliced data
                    sequence = str(record['sequence'][21:])  # cut off the index + 12 random nucleotide (leaving 130 nt)
                    quality = str(record['quality'][21:])  # cut off the same amount of quality scores

                    # read 2 sliced data
                    sequence2 = str(record2['sequence'][:142])  # cut off the 8 random nucleotide sequence at end
                    quality2 = str(record2['quality'][:142])  # cut off the same amount of quality scores

                else:  # is the in-house AAFC barcodes.
                    index_name, slicer = parse_index(record, barcodes)
                    barcodes[index_name][1] += 1

                    # read 1 sliced data
                    sequence = str(record['sequence'][slicer:])  # cut off the index, leaving just TCGA at start
                    quality = str(record['quality'][slicer:])  # cut off the same amount of quality scores

                    # read 2 sliced data
                    sequence2 = str(record2['sequence'])  # no index on read2 so no need to slice sequence like read1
                    quality2 = str(record2['quality'])

                # double check that the run_ids match between read1 and read 2
                run_id = record['name'].split(' ')[0]
                run_id2 = record2['name'].split(' ')[0]
                read_id1 = record['name'].split(' ')[1]
                read_id2 = record2['name'].split(' ')[1]

                if run_id == run_id2:

                    if ':0:0' in run_id:  # In house data had an issue where coordinates were 0:0 for all lines.
                        x_coord = str(rand_x)
                        y_coord = str(rand_y)

                        # apply the random incremental integers to the coordinates instead of 0:0
                        run_id = run_id.replace(':0:0', f'{x_coord}:{y_coord}')
                        record['name'] = f'{run_id} {read_id1}'  # recreate the record name for R1 with coordinates
                        record2['name'] = f'{run_id} {read_id2}'  # recreate the record name for R2 with coordinates

                        rand_x += 1  # increment so next coordinate isn't the same.
                        rand_y += 1  # increment so next coordinate isn't the same.

                    # add the record information to the samples dictionary as a string
                    fastq1 = f"{record['name']}\n{sequence}\n+\n{quality}\n"
                    fastq2 = f"{record2['name']}\n{sequence2}\n+\n{quality2}\n"

                    # sometimes they match to index not actually in the sample sheet, so weed those out
                    if index_name in samples:
                        # add them to the dictionaries value[2] (the list)
                        samples[index_name][2].append(fastq1)
                        samples2[index_name][2].append(fastq2)

                    else:  # assign it to unmatched instead
                        samples['unmatched'][2].append(fastq1)
                        samples2['unmatched'][2].append(fastq2)

                    # reset your lines for next sequence record in fastq
                    lines1 = []
                    lines2 = []

                else:
                    print('your fastq file formats are not matching.  Could be issues with fastq files')
                    continue

    # write the files for the sample dictionary entries for read 1
    fq_write_lock(lock, read1=samples, read2=samples2)

    return barcodes


def write_outstats(st_time, fq_path, bcode_file, bcode_dict, sample_dict):
    """
    :param st_time:  time that the program started
    :param fq_path:  The -f argument call, which includes the path in which to create and save the log
    :param bcode_file:  the -b argument call.  Add to the header info on what was run
    :param bcode_dict:  Barcode dictionary made from barcode file, includes barcode and counter as values for keys
    :param sample_dict:  sample dictionary made from samplesheet, used for adding number of bcodes found per sample
    :return:  writes/appends the command, start and finish, and final output of demultiplex script to log file
    """
    #  write total demultiplex stats for the script (no split into reference genomes)
    write_name = f'{os.path.split(fq_path)[0]}/demultiplex.output.log'

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

                    # add the total number from bcode dict to the sample dict, only if that barcode is in the sample dictionary
                    if bcode_key in sample_dict:
                        sample_dict[bcode_key][2].append(bcode_val[1])
                    else:
                        # the barcode found doesnt match to a sample in dict, don't add it
                        continue

            percent_partitioned = (total_pcount / total_rcount) * 100
            reads = f'\nTotal reads in PE files: {total_rcount}\n'
            partitioned = f'Total partitioned reads: {total_pcount} \n'
            percent_p = f'Percentage partitioned : {percent_partitioned:.2f}'
            outfile.write(reads)
            outfile.write(partitioned)
            outfile.write(percent_p)

        # write the outstats of each demultiplexed reference folder. Get just the directories from sample_dict[index][0]
        ref_dirs = []
        for sam_key, sam_val in sample_dict.items():
            if sam_val[0] not in ref_dirs:
                ref_dirs.append(sam_val[0])

        # write the stats for each split reference directories log
        for r_dir in ref_dirs:
            if r_dir != 'unmatched':
                ref_dir_path = f"{os.path.split(fq_path)[0]}/{r_dir}/log/{r_dir}.demultiplex.stats"
                header_line = "Sample_name\tIndex_name\ttotal_entries\n"

                with open(ref_dir_path, 'w') as ref_out:
                    ref_out.write(header_line)

                    for sam_key, sam_val in sample_dict.items():
                        if sam_val[0] == r_dir:
                            log_line = f"{sam_val[1]}\t{sam_key}\t{sam_val[2][0]}\n"
                            ref_out.write(log_line)


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
    write_outstats(st_time=start_time, fq_path=fastq_file1, bcode_file=barcode_file, bcode_dict=None, sample_dict=None)

    # check to see if demultiplex folder exists for this sample, if not, create it
    path_name, _ = os.path.split(fastq_file1)

    # create the demultiplex stats dictionary, and generate the demultiplex and log folders (if not already created)
    demultiplexed_stats, sample_paths = make_dicts(barcode_file, sample_sheet)

    for sample_key, sample_val in sample_paths.items():
        write_path = f"{path_name}/{sample_val[0]}/demultiplex/"
        Path(write_path).mkdir(parents=True, exist_ok=True)

        if sample_val[0] != 'unmatched':
            log_path = f"{path_name}/{sample_val[0]}/log"
            Path(log_path).mkdir(parents=True, exist_ok=True)

    # get a list of the chunked data to demultiplex based off fastq_R1 path/name
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

    ####################################  write out final stats  ###############################################

    # write the stdout for the barcodes founds in the file.  total reads, etc.
    write_outstats(st_time=None, fq_path=fastq_file1, bcode_file=barcode_file, bcode_dict=demultiplexed_stats,
                   sample_dict=sample_paths)
