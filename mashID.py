import os
from argparse import ArgumentParser
from mashID_methods import Methods
import pkg_resources
import shutil
from multiprocessing import cpu_count
from psutil import virtual_memory


__author__ = 'duceppemo'
__version__ = '0.1'


class Identification(object):
    def __init__(self, args):
        # Command line arguments

        # Performance
        self.cpu = args.threads
        self.parallel = args.parallel
        self.mem = args.memory

        # I/O
        self.input = os.path.abspath(args.input)
        self.output_folder = os.path.abspath(args.output)
        if args.database:
            self.mash_db = os.path.abspath(args.databse)
        else:  # use the default Mycobacteria DB
            self.mash_db = pkg_resources.resource_filename('dependencies', 'mycobacteria_mash_sketches.msh')

        # Data
        self.sample_dict = dict()

        # Run
        self.run()

    def run(self):
        # Create output folders
        Methods.make_folder(self.output_folder)

        # Get input files and place info in dictionary
        print('Getting input file(s) stats...')
        if os.path.isdir(self.input):  # single or multiple files inside folder
            self.sample_dict = Methods.get_files(self.input)
        else:  # is a single file
            sample = os.path.basename(self.input).split('.')[0].replace('_pass', '')
            self.sample_dict = {sample: [os.path.realpath(self.input)]}  # Follow symbolic links

        # Merge paired-end fastq files is present and update sample dict
        tmp_folder = self.output_folder + '/tmp'
        Methods.make_folder(tmp_folder)
        self.sample_dict = Methods.merge_fastq_pair(self.sample_dict, tmp_folder)

        # Get file stats (number of reads/contigs and bp)
        Methods.get_stats_parallel(self.sample_dict, self.mem, self.cpu, self.parallel)

        # Screen samples and create summary report
        print('Identifying samples...')
        samples_df = Methods.mash_screen_parallel(self.mash_db, self.output_folder, self.sample_dict)

        # Print summary report to terminal
        print('\nIdentification results:\n')
        print(samples_df.to_string(index=False, justify='left'))

        # Write output file
        if len(self.sample_dict) > 1:
            output_tsv = self.output_folder + '/summary_mashID.tsv'
            samples_df.to_csv(output_tsv, sep="\t", index=False)

        # Delete tmp folder (needed to merged paired-end fastq files, if present)
        shutil.rmtree(tmp_folder)


if __name__ == "__main__":
    max_cpu = cpu_count()
    max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB

    parser = ArgumentParser(description='Species identification from NGS data using Mash.')
    parser.add_argument('-i', '--input', metavar='/input/folder/',
                        help='Input directory with fastq/fasta files or a single fastq/fasta file, gzipped or not.',
                        type=str, required=True)
    parser.add_argument('-o', '--output', metavar='/output/folder/',
                        help='Output directory',
                        type=str, required=True)
    parser.add_argument('-d', '--database', metavar='/path/to/mash_databse.msh',
                        help='Mash sketch database. Will run a Mycobacteria mash database by default.',
                        type=str, required=False)
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False,
                        type=int, default=max_cpu,
                        help='Number of threads. Default is maximum available({}). Optional.'.format(max_cpu))
    parser.add_argument('-p', '--parallel', metavar='2',
                        required=False,
                        type=int, default=2,
                        help='Number of samples to process in parallel. Default is 2. Optional.')
    parser.add_argument('-m', '--memory', metavar=str(max_mem),
                        required=False,
                        type=int, default=max_mem,
                        help='Memory in GB. Default is 85%% of total memory ({})'.format(max_mem))
    parser.add_argument('-v', '--version', action='version',
                        version=f'{os.path.basename(__file__)}: version {__version__}')

    arguments = parser.parse_args()
    Identification(arguments)
