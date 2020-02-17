#!usr/bin/env python3
"""
This script is used to demultiplex paired end fastq format GBS files into the seperate GBS index
files based on the in-line barcodes annealed to the DNA fragments in Read1 files only
(The GBS data to be demuxed is produced by the iGenomics Riptide Kit, and sequenced 2x150 PE on Illumina platform)
To run this script individually and not part of a pipeline use the following example command:
$ python3 scripts/riptide.py -f path/to/read_1.fastq -s path/to/samplesheet.txt -b path/to/barcode_file.txt
Author:  Brian James
email:  brian.james4@canada.ca
"""

import re
import argparse
import os
from datetime import datetime
import sqlite3
import gzip
import shutil


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
            sample_dict[sam[1]] = sam[2]
        sample_dict['unmatched'] = 'unmatched'

    return barcode_dict, sample_dict


def process(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}


def db_connection(db_path):
    """
    create a temp database based on the path above, and return a connection after creating some tables in the database
    :param db_path:
    :return: database connection (db.cursor())
    """
    # create a temp SQLite database connection that houses the temp data for writing new fastq files
    sql = sqlite3.connect(f"{db_path}/temp.db")
    sql.isolation_level = None  # turn on auto commit
    connection = sql.cursor()  # create connection to database
    connection.execute("PRAGMA cache_size=100000")  # add in parameters for max batch sizes of 100,000

    return connection


def create_tempsql(connection, barcode_dict):
    """create the temp database"""

    conn.execute('''CREATE TABLE IF NOT EXISTS records (
        id INTEGER PRIMARY KEY AUTOINCREMENT, 
        record_id TEXT UNIQUE, 
        index_name TEXT)''')

    # create 192 tables for each index in the barcode list, with room for each read
    for ind_name, ind_bcode in barcode_dict.items():
        table_make = f'CREATE TABLE IF NOT EXISTS {ind_name} (record_id1 TEXT UNIQUE, record_val1 TEXT, sequence1 TEXT, ' \
                     f'quality1 TEXT, record_val2 TEXT, sequence2 TEXT, quality2 TEXT)'

        connection.execute(table_make)


def unzipper(fastq):
    """
    read in a gzipped file, and return the unzipped version names
    :param fastq:  the file from the arguments in gzip format (read 1 only)
    """
    # split up the file name and extensions
    fastq_2 = fastq.replace("R1", "R2")
    path_name, base_name = os.path.split(fastq)
    file_handle = os.path.splitext(base_name)[0]

    unzipped = f'{path_name}/{file_handle}'
    unzipped2 = unzipped.replace("R1", "R2")

    with gzip.open(fastq, 'rb') as f_in:
        with open(unzipped, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    os.remove(fastq)

    with gzip.open(fastq_2, 'rb') as f2_in:
        with open(unzipped2, 'wb') as f2_out:
            shutil.copyfileobj(f2_in, f2_out)

    os.remove(fastq_2)

    return unzipped, unzipped2


def parse_read1(seq_record, bcode_dict):
    """
    :param seq_record: the record parsed from SeqIO.parse, houses info from fastq file
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

        # only keep the matches that are the same length as the index (if multiple hits, reduce to 1)
        for key, val in gbs_index.items():
            return_name.append(key)

    if len(return_name) == 0:
        return 'unmatched'
    else:
        return return_name[0]


def insert_fqdata(fastq, db_path, read_num, barcode_dict):
    """
    read through the fastq files and insert the data into the database, done in batches of 100,000 lines at a time.
    :param fastq:  the fastq file you would like to read and pass into database
    :param db: database connection
    :param read_num: value (1, 2):  1=fastq_R1  2=fastq_R2
    :param barcode_dict: the dictionary for the barcodes in the samplesheet
    :return: updated barcode dict
    """
    db = db_connection(db_path)
    db.execute("BEGIN TRANSACTION")  # start transaction so we only insert data when we reach size size

    with open(fastq, 'r') as fq:
        lines = []

        for line in fq:
            lines.append(line.rstrip())
            if len(lines) == 4:
                record = process(lines)

                if read_num == 1:
                    index_name = parse_read1(record, barcode_dict)
                    barcode_dict[index_name][1] += 1  # add instance of the index being counted into the barcode dict

                    #sequence = record['sequence'][20:]  # cut off the index, leaving just TCGA at start of sequence
                    sequence = record['sequence']
                    #quality = record['quality'][20:]  # cut off the same amount of quality scores
                    quality = record['quality']
                    run_id = str(record['name']).split(' ')[0]
                    run_val = str(record['name']).split(' ')[1]

                    # insert statements for the read 1, inserting data into read1 table and unmatched table
                    insert_data = f'''INSERT INTO {index_name} (record_id1, record_val1, sequence1, quality1) ''' \
                                  f'''VALUES ('{run_id}', '{run_val}', '{sequence}', '{quality}')'''
                    unmatched_data = f'''INSERT INTO unmatched (record_id1, record_val1, sequence1, quality1) ''' \
                                     f'''VALUES ('{run_id}', '{run_val}', '{record["sequence"]}', '{record["quality"]}')'''
                    insert_index = f'INSERT INTO records (record_id, index_name) VALUES ("{run_id}", "{index_name}")'

                    if index_name != 'unmatched':  # index was found, insert into proper table
                        db.execute(insert_data)
                        db.execute(insert_index)

                    else:  # add the data to the unmatched table
                        db.execute(unmatched_data)
                        db.execute(insert_index)

                else:  # its read 2, so no index_names are available
                    run_id = str(record['name']).split(' ')[0]  # split the record run_id to help match
                    run_val = str(record['name']).split(' ')[1]
                    sequence = str(record['sequence'])
                    quality = str(record['quality'])
                    query = f'''SELECT index_name FROM records WHERE record_id="{run_id}"'''
                    fetched = conn.execute(query).fetchall()
                    index_name = fetched[0][0]

                    insert_data = f'UPDATE {index_name} SET record_val2="{run_val}",' \
                                  f' sequence2="{sequence}", quality2="{quality}" WHERE record_id1="{run_id}" '
                    db.execute(insert_data)

                lines = []  # reset the lines at end of loop

    db.execute('''END TRANSACTION''')
    db.close()


def savefile_handle(fastq_path, sample_name):
    """
    :param fastq_path:  sample1_R1.fastq, sample1_R2.fastq
    :param sample_name: pass in sample name (from dictionary) that will match index to sample name.
    :return:  SK-GBD-000001.1.fastq, SK-GBD-000001.2.fastq, SK-GBD-000002.1.fastq, SK-GBD-000002.2.fastq, etc.
    """
    # split up the fastq file extension from name, add the index, then return the save file name
    path_name, base_name = os.path.split(fastq_path)
    extension = os.path.splitext(base_name)[1]
    file_handle = os.path.splitext(base_name)[0]

    if 'R1' in file_handle:
        save_handle = f'{path_name}/demultiplex/{sample_name}.1{extension}'
    else:
        save_handle = f'{path_name}/demultiplex/{sample_name}.2{extension}'

    # check to see if demultiplex folder exists for this sample, if not, create it
    if os.path.isdir(f'{path_name}/demultiplex/'):
        pass
    else:
        os.mkdir(f'{path_name}/demultiplex/')

    return save_handle


def write_fq(save_handle, db_path, type_fq, index_name=None):
    """write the query out to a seperate .fastq file for each index in the file
    :param save_handle: the name of the file you want to write out
    :param db: database connection
    :param type_fq: values ('matched', 'unmatched1', 'unmatched2')
    :param index_name:  the index to query on in database
    :return writes out the fastq file into the demultiplexed folder
    """
    db = db_connection(db_path)
    # combine these into 1 query, write open both files, and write both files line by line (not with open())
    if type_fq == 'unmatched':
        save_handle2 = save_handle.replace(".1.", ".2.")
        fq_data = list(db.execute('''SELECT * FROM unmatched''').fetchall())

        r1 = open(save_handle, "w")
        r2 = open(save_handle2, "w")

        for entry in fq_data:
            write_data1 = [f'{entry[0]} {entry[1]}\n', f'{entry[2]}\n', f'+\n', f'{entry[3]}\n']
            write_data2 = [f'{entry[0]} {entry[4]}\n', f'{entry[5]}\n', f'+\n', f'{entry[6]}\n']
            r1.writelines(write_data1)
            r2.writelines(write_data2)

        r1.close()
        r2.close()

    # unmatched fastq files
    else:
        save_handle2 = save_handle.replace(".1.", ".2.")
        sample_query = f'''SELECT * FROM {index_name}'''
        fq_data = list(db.execute(sample_query).fetchall())

        r1 = open(save_handle, "w")
        r2 = open(save_handle2, "w")

        for entry in fq_data:
            write_data1 = [f'{entry[0]} {entry[1]}\n', f'{entry[2]}\n', f'+\n', f'{entry[3]}\n']
            write_data2 = [f'{entry[0]} {entry[4]}\n', f'{entry[5]}\n', f'+\n', f'{entry[6]}\n']
            r1.writelines(write_data1)
            r2.writelines(write_data2)

        r1.close()
        r2.close()
    db.close()

def write_outstats(st_time, fastq_path, bcode_file, bcode_dict):
    """
    :param st_time:  time that the program started
    :param fastq_path:  The -f argument call, which includes the path in which to create and save the log
    :param bcode_file:  the -b argument call.  Add to the header info on what was run
    :param bcode_dict:  Barcode dictionary made from barcode file, includes barcode and counter as values for keys
    :return:  writes/appends the command, start and finish, and final output of demultiplex script to log file
    """

    if os.path.isdir(f'{os.path.split(fastq_path)[0]}/log'):
        pass
    else:
        os.mkdir(f'{os.path.split(fastq_path)[0]}/log')

    write_name = f'{os.path.split(fastq_path)[0]}/log/demultiplex.output.log'

    command = f'Working Directory:  {os.getcwd()}\n'\
              f'Running Command:\n\t$ python3 scripts/fastq_demultiplex.py -f {fastq_path} -b {bcode_file}\n'

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
    fastq_file2 = fastq_file1.replace("R1", "R2")   # read 2.fastq file name
    barcode_file = args.barcode                     # barcode file (tab delimited, 2 headers)
    sample_file = args.sample                       # sample config file (tab delimited, 3 headers)
    start_time = datetime.now()                     # get the start time for the program starting

    # determine if the file being read is in gzip format
    if '.gz' in fastq_file1:
        fastq_file1, fastq_file2 = unzipper(fastq_file1)

    # write header and info for stats output/log file
    write_outstats(st_time=start_time, fastq_path=fastq_file1, bcode_file=barcode_file, bcode_dict=None)

    # get barcode and sample file into a dictionary format for passing into parse_barcode() function
    barcodes, samples = make_dicts(barcode_file, sample_file)

    # create a temp SQLite database connection that houses the temp data for writing new fastq files
    db_path = os.path.split(fastq_file1)[0]
    conn = db_connection(db_path)
    create_tempsql(conn, barcodes)
    conn.close()

    #################################### demultiplex read1 and read2 fastq data ########################################

    insert_fqdata(fastq_file1, db_path, read_num=1, barcode_dict=barcodes)  # insert the read1 data into sqlite db

    insert_fqdata(fastq_file2, db_path, read_num=2, barcode_dict=barcodes)  # insert the read2 data into sqlite db

    # parse the database and split up the data into separate fastq files based on index_name
    for key, val in samples.items():
        # write read1 and 2 per index in the sample sheet
        save_name = savefile_handle(fastq_file1, sample_name=val)
        write_fq(save_name, db_path, type_fq='matched', index_name=key)

    # generate the name for unmatched files an write out the unmatched data
    unmatched_data = savefile_handle(fastq_file1, sample_name='unmatched')
    write_fq(unmatched_data, db_path, type_fq='unmatched')

    ####################################  write out stats, close program ###############################################

    # write the stdout for the barcodes founds in the file.  total reads, etc.
    write_outstats(st_time=None, fastq_path=fastq_file1, bcode_file=barcode_file, bcode_dict=barcodes)
    conn.close()

    os.remove(f"{os.path.split(fastq_file1)[0]}/temp.db")
