import subprocess
import os
import pathlib
import shutil
from io import StringIO
from io import TextIOWrapper
import pandas as pd
from concurrent import futures
import gzip
from itertools import groupby
from multiprocessing import cpu_count
from psutil import virtual_memory
import sys


class Methods(object):
    accepted_extensions = ['.fq', '.fq.gz',
                           '.fastq', '.fastq.gz',
                           '.fasta', '.fasta.gz',
                           '.fa', '.fa.gz',
                           '.fna', '.fna.gz']

    @staticmethod
    def check_cpus(requested_cpu, n_proc):
        total_cpu = cpu_count()

        if 1 > requested_cpu > total_cpu:
            requested_cpu = total_cpu
            sys.stderr.write("Number of threads was set to {}".format(requested_cpu))
        if 1 > n_proc > total_cpu:
            n_proc = total_cpu
            sys.stderr.write("Number of samples to parallel process was set to {}".format(total_cpu))

        return requested_cpu, n_proc

    @staticmethod
    def check_mem(requested_mem):
        max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB
        if requested_mem:
            if requested_mem > max_mem:
                requested_mem = max_mem
                sys.stderr.write("Requested memory was set higher than available system memory ({})".format(max_mem))
                sys.stderr.write("Memory was set to {}".format(requested_mem))
        else:
            requested_mem = max_mem

        return requested_mem

    @staticmethod
    def check_identity_range(my_range):
        if 0 < my_range > 1:
            raise Exception('identity value must be between 0 and 1.')

    @staticmethod
    def check_p_value(my_p_value):
        if 0 < my_p_value > 1:
            raise Exception('P-value must be between 0 and 1.')

    @staticmethod
    def check_input(my_input):
        if not os.path.isdir(my_input):
            raise Exception('Please provide a folder containing fastq files as input.')
        if not os.path.exists(my_input):
            raise Exception('The provided input folder does not exist.')

    @staticmethod
    def make_folder(folder):
        # Will create parent directories if they don't exist and will not return error if already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def get_files(my_input):
        sample_dict = dict()

        if os.path.isdir(my_input):  # Input is a folder
            # Look for input sequence files recursively
            for root, directories, filenames in os.walk(my_input):
                for filename in filenames:
                    if filename.endswith(tuple(Methods.accepted_extensions)):  # accept a tuple or string
                        file_path = os.path.join(root, filename)
                        file_path = os.path.realpath(file_path)  # follow symbolic links

                        # Get sample name
                        # sample = filename.split('.')[0].replace('_pass', '').replace('_filtered', '')
                        sample = filename.split('.')[0].split('_')[0]
                        if filename.endswith('.gz'):
                            sample = sample.split('.')[0]

                        # Get total reads and bp for fastq/fasta

                        # Create dictionary entries
                        if sample not in sample_dict.keys():
                            sample_dict[sample] = dict()
                            sample_dict[sample]['path'] = list()
                            sample_dict[sample]['reads'] = list()
                            sample_dict[sample]['bp'] = list()

                        sample_dict[sample]['path'].append(file_path)
        elif os.path.isfile(my_input):  # Input is a file
            sample = os.path.basename(my_input).split('.')[0].replace('_pass', '')
            sample_dict[sample] = dict()
            sample_dict[sample] = {'path': [os.path.realpath(my_input)]}  # Follow symbolic links
        else:
            raise Exception('Hmmm... something went terribly wrong!')

        if not sample_dict:
            raise Exception('Sample dictionary empty!')

        return sample_dict

    @staticmethod
    def merge_fastq_pair(sample_dict, output_folder):
        # Check if multiple files per sample (if paired-end)
        for sample, info_dict in sample_dict.items():
            if len(info_dict['path']) == 2:  # Paired-end
                # Figure out if input file is gzipped or not
                file_ext = info_dict['path'][0].split('.')[-1]
                if file_ext == 'gz':
                    file_ext = info_dict['path'][0].split('.')[-2]
                    merged_file = output_folder + '/' + sample + '.' + file_ext + '.gz'
                else:
                    merged_file = output_folder + '/' + sample + '.' + file_ext

                with open(merged_file, 'wb') as wfp:
                    for seq_file in info_dict['path']:
                        with open(seq_file, 'rb') as rfp:
                            shutil.copyfileobj(rfp, wfp)
            else:
                try:
                    # just create symbolic link
                    os.symlink(info_dict['path'][0], output_folder +
                               '/' + os.path.basename(info_dict['path'][0]))
                except FileExistsError:
                    pass

        # Update sample dict with new merged files
        return Methods.get_files(output_folder)

    @staticmethod
    def mash_screen(sample, mash_db, sample_path, output_folder, identity, p_value, cpu, n_hit):
        print('\t{}'.format(sample))

        cmd = ["mash", "screen",
               '-i', str(identity),
               '-v', str(p_value),
               '-w',  # to reduce redundancy in output
               '-p', str(cpu),
               mash_db,
               sample_path]

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

        my_stringio = StringIO(p.communicate()[0].decode('utf-8'))
        # 'Identity', 'Shared-Hashes', 'Median-Multiplicity', 'P-Value', 'Query-ID', 'Query-Comment'
        df = pd.read_csv(my_stringio, sep="\t", names=['Identity', 'Shared-Hashes', 'Median-Multiplicity',
                                                       'P-Value', 'Query-ID', 'Query-Comment'])
        df = df.sort_values(by=['Identity'], kind='mergesort', ascending=False)
        df = df.head(n=n_hit)  # Only keep the top 5

        df = df.reset_index(drop=True)  # Renumber indexes
        df.index = df.index + 1  # Change Index column to starts at 1 instead of 0
        df.index.name = 'Rank'  # Change index column name

        # Reformat Query-ID to only keep the accession number
        df['Query-ID'] = list(map(lambda x: '.'.join(os.path.basename(x).split(".")[:-1]), df['Query-ID'].tolist()))

        # Write output file
        output_tsv = output_folder + '/' + sample + "_mashID.tsv"
        df.to_csv(output_tsv, sep="\t")

        return sample, df

    @staticmethod
    def mash_screen_parallel(mash_db, output_folder, sample_dict, identity, p_value, cpu, parallel, n_hit):
        df = pd.DataFrame(columns=['Sample', 'Reads', 'Length', 'Identity', 'Shared-Hashes',
                                   'Median-Multiplicity', 'P-Value', 'Query-ID'])

        with futures.ThreadPoolExecutor(max_workers=parallel) as executor:
            args = ((sample, mash_db, info_dict['path'][0], output_folder,
                     identity, p_value, int(cpu / parallel), n_hit) for sample, info_dict in sample_dict.items())
            for sample, my_df in executor.map(lambda x: Methods.mash_screen(*x), args):
                # Only keep best hit
                my_df = my_df.head(n=1)

                try:
                    ident = my_df['Identity'].iloc[0]
                    hashes = my_df['Shared-Hashes'].iloc[0]
                    mult = my_df['Median-Multiplicity'].iloc[0]
                    p_value = my_df['P-Value'].iloc[0]
                    query_id = my_df['Query-ID'].iloc[0]
                    comment = my_df['Query-Comment'].iloc[0]  # Get the ID
                    org_id = Methods.species_from_header(comment)  # Change name
                except IndexError:
                    ident = hashes = mult = p_value = query_id = 'NA'
                    org_id = 'No significant hit in database'

                results_dict = {'Sample': sample,
                                'Reads': sample_dict[sample]['reads'][0],
                                'Length': sample_dict[sample]['bp'][0],
                                'Identity': [ident],
                                'Shared-Hashes': [hashes],
                                'Median-Multiplicity': [mult],
                                'P-Value': [p_value],
                                'Query-ID': [query_id],
                                'Query-Comment': [org_id]
                                }
                new_df = pd.DataFrame.from_dict(results_dict)
                df = pd.concat([df, pd.DataFrame.from_dict(results_dict)], axis='index', ignore_index=True)
        # Sort df by sample
        df.sort_values(by=['Sample'], axis='index', ascending=True, inplace=True, ignore_index=True)

        # rename some columns
        df.rename(columns={'Identity': '%-Identity', 'Query-ID': 'Accession', 'Query-Comment': 'Identification'},
                  inplace=True)

        return df

    @staticmethod
    def species_from_header(header):
        if header.startswith('['):
            header_list = header.split()[2:]
        else:
            header_list = header.split()
        genus = header_list[1]
        species = header_list[2]

        if 'variant' in header:
            variant = header_list[4]
            new_name = '{} {} variant {}'.format(genus, species, variant)
        elif any(x in header for x in ['subsp', 'sub']):
            subsp = header_list[4]
            new_name = '{} {} subsp. {}'.format(genus, species, subsp)
        else:
            new_name = '{} {}'.format(genus, species)

        return new_name

    @staticmethod
    def get_read_bp(seq_file):
        total_reads = 0
        total_bp = 0
        with gzip.open(seq_file, 'rb', 1024 * 1024) if seq_file.endswith('.gz') \
                else open(seq_file, 'r', 1024 * 1024) as f:
            if any(x in seq_file for x in ['.fq', '.fastq']):  # if fastq
                counter = 0
                for line in f:
                    counter += 1
                    if counter == 2:
                        line = line.rstrip()
                        total_reads += 1
                        total_bp += len(line)
                    if counter == 4:
                        counter = 0
            else:  # if fasta
                # Need to be tested
                # https://www.biostars.org/p/710/#120760
                faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
                total_reads = len([x for x in faiter])
                total_bp = sum([len("".join(s.strip() for s in faiter.__next__())) for header in faiter])

        return total_reads, total_bp

    @staticmethod
    def get_stats(seq_file, sample, mem, cpu):
        cmd = ['stats.sh',
               '-Xmx{}g'.format(mem),
               'in={}'.format(seq_file),
               'threads={}'.format(cpu),
               'format=2']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

        reads = 0
        bp = 0
        for line in TextIOWrapper(p.stdout, encoding="utf-8"):
            line = line.strip()
            if 'Main genome contig total' in line:
                reads = line.split(':')[1].split('\t')[1].strip()
            if 'Main genome contig sequence total' in line:
                bp = line.split(':')[1].split('\t')[1].split()[0]
        return sample, reads, bp

    @staticmethod
    def get_stats_parallel(sample_dict, mem, cpu, parallel):
        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((info_dict['path'][0], sample, int(mem / parallel), int(cpu / parallel))
                    for sample, info_dict in sample_dict.items())
            for sample, reads, bp in executor.map(lambda x: Methods.get_stats(*x), args):
                sample_dict[sample]['reads'].append(reads)
                sample_dict[sample]['bp'].append(bp)
