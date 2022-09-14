import os
from argparse import ArgumentParser
from mashID_methods import Methods
import pkg_resources
import shutil


__author__ = 'duceppemo'
__version__ = '0.1'


class Identification(object):
    def __init__(self, args):
        # Command line arguments
        self.input = os.path.abspath(args.input)
        self.output_folder = os.path.abspath(args.output)
        if args.mash_db:
            self.mash_db = os.path.abspath(args.mash_db)
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
        if os.path.isdir(self.input):
            self.sample_dict = Methods.get_files(self.input)
        else:  # is a file
            sample = os.path.basename(self.input).split('.')[0].replace('_pass', '')
            self.sample_dict = {sample: [os.path.realpath(self.input)]}  # Follow symbolic links

        # Merge paired-end fastq files is present and update sample dict
        tmp_folder = self.output_folder + '/tmp'
        Methods.make_folder(tmp_folder)
        self.sample_dict = Methods.merge_fastq_pair(self.sample_dict, tmp_folder)

        # Screen samples
        print('Identifying samples...')
        samples_df = Methods.mash_screen_parallel(self.mash_db, self.output_folder, self.sample_dict)
        print('\nIdentification results:\n')
        print(samples_df.to_string(index=False, justify='left'))

        # Write output file
        if len(self.sample_dict) > 1:
            output_tsv = self.output_folder + '/summary_mashID.tsv'
            samples_df.to_csv(output_tsv, sep="\t", index=False)

        # Delete tmp folder
        shutil.rmtree(tmp_folder)


if __name__ == "__main__":
    parser = ArgumentParser(description='Species identification from NGS data using Mash.')
    parser.add_argument('-i', '--input', metavar='/input/folder/',
                        help='Input directory with fastq/fasta files or a single fastq/fasta file, gzipped or not.',
                        type=str, required=True)
    parser.add_argument('-o', '--output', metavar='/output/folder/',
                        help='Output directory',
                        type=str, required=True)
    parser.add_argument('-m', '--mash-db', metavar='/path/to/file.msh',
                        help='Mash sketch database. Will run a Mycobacteria mash database by default.',
                        type=str, required=False)
    parser.add_argument('-v', '--version', action='version',
                        version=f'{os.path.basename(__file__)}: version {__version__}')

    arguments = parser.parse_args()
    Identification(arguments)
