import subprocess
import os
import pathlib
import shutil
from io import StringIO
import pandas as pd
from concurrent import futures
from multiprocessing import cpu_count


class Methods(object):
    accepted_extensions = ['.fq', '.fq.gz',
                           '.fastq', '.fastq.gz',
                           '.fasta', '.fasta.gz',
                           '.fa', '.fa.gz',
                           '.fna', '.fna.gz']

    @staticmethod
    def make_folder(folder):
        # Will create parent directories if they don't exist and will not return error if already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def get_files(in_folder):
        sample_dict = dict()

        # Look for input sequence files recursively
        for root, directories, filenames in os.walk(in_folder):
            for filename in filenames:
                if filename.endswith(tuple(Methods.accepted_extensions)):  # accept a tuple or string
                    file_path = os.path.join(root, filename)
                    file_path = os.path.realpath(file_path)  # follow symbolic links
                    sample = filename.split('.')[0].replace('_pass', '').replace('_filtered', '')
                    if filename.endswith('.gz'):
                        sample = sample.split('.')[0]
                    if sample not in sample_dict:
                        sample_dict[sample] = []
                    sample_dict[sample].append(file_path)
        if not sample_dict:
            raise Exception('Sample dictionary empty!')

        return sample_dict

    @staticmethod
    def merge_fastq_pair(sample_dict, output_folder):
        # Check if multiple files per sample (if paired-end)
        for sample, seq_list in sample_dict.items():
            if len(seq_list) == 2:  # Paired-end
                # Figure out if input file is gzipped or not
                file_ext = seq_list[0].split('.')[-1]
                if file_ext == '.gz':
                    file_ext = seq_list[0].split('.')[-2]
                    merged_file = output_folder + '/' + sample + file_ext + '.gz'
                else:
                    merged_file = output_folder + '/' + sample + file_ext

                with open(merged_file, 'wb') as wfp:
                    for seq_file in seq_list:
                        with open(seq_file, 'rb') as rfp:
                            shutil.copyfileobj(rfp, wfp)
            else:
                try:
                    # just create symbolic link
                    os.symlink(seq_list[0], output_folder + '/' + os.path.basename(seq_list[0]))
                except FileExistsError:
                    pass

        # Update sample dict with new merged files
        return Methods.get_files(output_folder)

    @staticmethod
    def mash_screen(sample, mash_db, sample_path, output_folder):
        print('\t{}'.format(sample))

        cmd = ["mash", "screen",
               '-i', str(0.9),
               '-v', str(0.05),
               '-w',  # to reduce redundancy in output
               mash_db,
               sample_path]

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

        my_stringio = StringIO(p.communicate()[0].decode('utf-8'))
        # 'Identity', 'Shared-Hashes', 'Median-Multiplicity', 'P-Value', 'Query-ID', 'Query-Comment'
        df = pd.read_csv(my_stringio, sep="\t", names=['Identity', 'Shared-Hashes', 'Median-Multiplicity',
                                                       'P-Value', 'Query-ID', 'Query-Comment'])
        df = df.sort_values(by=['Identity'], kind='mergesort', ascending=False)
        df = df.head(n=5)  # Only keep the top 5

        df = df.reset_index(drop=True)  # Renumber indexes
        df.index = df.index + 1  # Change Index column to starts at 1 instead of 0
        df.index.name = 'Rank'  # Change index column name

        # Reformat Query-ID
        df['Query-ID'] = list(map(
            lambda x: '.'.join(os.path.basename(x).replace("_", ' ').split(".")[:-1]).split("strain")[0],
            df['Query-ID'].tolist()))

        # Drop last "Comment" column
        df = df.iloc[:, :-1]

        # Write output file
        output_tsv = output_folder + '/' + sample + "_mashID.tsv"
        df.to_csv(output_tsv, sep="\t")

        return sample, df

    @staticmethod
    def mash_screen_parallel(mash_db, output_folder, sample_dict):
        df = pd.DataFrame(columns=['Sample', 'Identity', 'Shared-Hashes', 'P-Value', 'Query-ID'])

        with futures.ThreadPoolExecutor(max_workers=cpu_count()) as executor:
            args = ((sample, mash_db, seq_list[0], output_folder)
                    for sample, seq_list in sample_dict.items())
            for sample, my_df in executor.map(lambda x: Methods.mash_screen(*x), args):
                my_df = my_df.head(n=1)
                results_dict = {'Sample': sample,
                                'Identity': my_df['Identity'],
                                'Shared-Hashes': my_df['Shared-Hashes'],
                                'P-Value': my_df['P-Value'],
                                'Query-ID': my_df['Query-ID']}
                df = pd.concat([df, pd.DataFrame.from_dict(results_dict)], axis='index', ignore_index=True)
        # Sort df by sample
        df.sort_values(by=['Sample'], axis='index', ascending=True, inplace=True, ignore_index=True)

        return df
