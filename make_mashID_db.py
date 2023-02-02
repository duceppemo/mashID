from argparse import ArgumentParser
import os
import pathlib
from multiprocessing import cpu_count
import subprocess


__author__ = 'duceppemo'
__version__ = '0.1.1'


class MakeDB(object):
    accepted_extensions = ['.fna', '.fna.gz',
                           '.fa', '.fa.gz',
                           '.fasta', '.fasta.gz']

    def __init__(self, args):
        # Check input
        self.input = args.input
        self.output = args.output
        self.prefix = args.prefix
        self.cpu = args.threads
        self.sketch_size = args.sketch_size
        self.kmer_size = args.kmer_size

        # Checks
        if self.cpu > cpu_count():
            self.cpu = cpu_count()
        elif self.cpu < 1:
            self.cpu = 1

        if not os.path.isdir(self.input):
            raise Exception('Please provide a folder as input.')

        if not os.path.exists(self.output):
            self.make_folder(self.output)

        if 1 < self.kmer_size > 32:
            raise Exception('Kmer size must be between 1-32, inclusively.')

        if self.sketch_size < 1:
            raise Exception('At least one sketch per sample must be used.')

        # Write a temporary file with all fasta
        fasta_list_file = self.output + '/sample_list.txt'
        self.list_fasta(self.input, self.accepted_extensions, fasta_list_file)

        self.run_mash_sketch(fasta_list_file, self.output, self.prefix, self.cpu, self.sketch_size, self.kmer_size)

    @staticmethod
    def list_fasta(input_folder, ext, list_file):
        input_list = list()

        # Look for input sequence files recursively
        for root, directories, filenames in os.walk(input_folder):
            for filename in filenames:
                if filename.endswith(tuple(ext)):  # accept a tuple or string
                    input_list.append(os.path.join(root, filename))

        if not input_list:  # If list file is empty
            raise Exception('No valid fasta files found in {}. Please make sure file ends with {}.'.format(
                input_folder, ','.join(ext)))

        with open(list_file, 'w') as f:
            for fasta_file in input_list:
                f.write(fasta_file + '\n')

    @staticmethod
    def make_folder(folder):
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def run_mash_sketch(input_list, output_folder, prefix, cpu, sketch_size, kmer_size):
        cmd = ['mash', 'sketch',
               '-p', str(cpu),
               '-l', input_list,
               '-o', output_folder + '/' + prefix,
               '-s', str(sketch_size),
               '-k', str(kmer_size)]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


if __name__ == "__main__":
    parser = ArgumentParser(description='Create a mash database.')
    parser.add_argument('-i', '--input', metavar='/input/folder/',
                        help='Input directory with fasta files, gzipped or not.',
                        type=str, required=True)
    parser.add_argument('-o', '--output', metavar='/output/folder/',
                        help='Output folder to save mash database. ".msh" will be added.',
                        type=str, required=True)
    parser.add_argument('-p', '--prefix', metavar='my_database',
                        help='Mash database name. ".msh" will be added.',
                        type=str, required=False, default='myDB')
    parser.add_argument('-t', '--threads', metavar='4',
                        help='Number of CPU. Default is 4.',
                        type=int, required=False, default=4)
    parser.add_argument('-s', '--sketch-size', metavar='1000',
                        help='Each sketch will have at most this many non-redundant min-hashes.',
                        type=int, required=False, default=1000)
    parser.add_argument('-k', '--kmer-size', metavar='21',
                        help='Hashes will be based on strings of this many nucleotides (1-32) . Default is 21.',
                        type=int, required=False, default=21)
    parser.add_argument('-v', '--version', action='version',
                        version=f'{os.path.basename(__file__)}: version {__version__}')

    # Get the arguments into an object
    arguments = parser.parse_args()

    MakeDB(arguments)
