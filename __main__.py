import os.path
import glob
import shutil

from trim import run_trim
from fastqc import run_fastqc, cut_rna_below_cutoff
from genome_align import align_to_genome, bin_rna_size
from unitas import run_unitas_annotation, merge_summary

from config import get_config_key, mkdir_if_not_exists, load_config

from argparse import ArgumentParser

def process_command(small_rna, adapter, front, anywhere, cutoff, quiet):
    print('==> Starting command Process')
    small_rna_path = small_rna

    if adapter is not None or front is not None or anywhere is not None:
        run_trim(small_rna_path, adapter, front, anywhere)
        small_rna_path = os.path.join(get_config_key('general', 'output_directory'), 'trimmed_rna.fastq')
    else:
        print('==> No adapter sequence provided, skipping trim step')

    if cutoff > 0:
        run_fastqc(small_rna_path, quiet=quiet)
        cut_rna_below_cutoff(small_rna_path, cutoff)
    else:
        # create the cut_sequences.fastq, even if fastqc wasn't run
        shutil.copy2(small_rna_path, os.path.join(get_config_key('general', 'output_directory'), 'cut_sequences.fastq'))

    print('==> Completed command Process')

def sort_command(genome, small_rna, min_length, max_length, quiet):
    print('==> Starting command Sort')
    new_fastq = align_to_genome(genome, small_rna, quiet=quiet)
    bin_rna_size(new_fastq, min_length, max_length)

    print('==> Completed command Sort')

def unitas_command(small_rna_path, species_name, ref_seqs, ref_seqs2, quiet):
    print('==> Starting command Unitas')
    UNITAS_OUTPUT = os.path.join(get_config_key('general', 'output_directory'), 'unitas')

    mkdir_if_not_exists(UNITAS_OUTPUT)

    for small_rna in glob.glob(os.path.join(small_rna_path, '*.fastq')):
        run_unitas_annotation(small_rna, species_name, ref_seqs, ref_seqs, quiet=quiet, unitas_output=UNITAS_OUTPUT)

    merge_summary()
    print('==> Completed command Unitas')

if __name__ == '__main__':
    parser = ArgumentParser(description='stuff')
    parser.add_argument('-q', '--quiet', help='Supress output from intermediate commands', action='store_true')
    parser.add_argument('-C', '--path-to-config', help='Path to the TOML format config file to use', default='config.toml')

    subparsers = parser.add_subparsers(dest='command')

    parser_process = subparsers.add_parser('process', help='Preprocessing for the RNA')
    parser_process.add_argument('-a', '--adapter', help='Sequence of the adapter to remove from the 3\' end')
    parser_process.add_argument('-g', '--front', help='Sequence of the adapter to remove from the 5\' end')
    parser_process.add_argument('-b', '--anywhere', help='Sequence of the adapters to remove from both ends')
    parser_process.add_argument('-c', '--cutoff', help='Quality cutoff to trin RNA sequences at', default=20, type=int)
    parser_process.add_argument('small_rna', help='Path to FASTQ containing the small RNA')

    parser_sort = subparsers.add_parser('sort', help='Find RNAs that align to a genome and sort them by length')
    parser_sort.add_argument('-l', '--min-length', help='Minimum length to bin', type=int)
    parser_sort.add_argument('-x', '--max-length', help='Maximum length to bin', type=int)
    parser_sort.add_argument('small_rna', help='Path to FASTQ containing the small RNA')
    parser_sort.add_argument('genome', help='Genome to align against')

    parser_unitas = subparsers.add_parser('unitas', help='Run unitas on split files and merge results')
    parser_unitas.add_argument('-r', '--refseq', help='References for use with unitas', nargs='*')
    parser_unitas.add_argument('-u', '--subsiquent-refseq', help='Refseqs to use in a second run', nargs='*')
    parser_unitas.add_argument('-s', '--species', help='Species to set in unitas arguments', default='x')
    parser_unitas.add_argument('path_to_rnas', help='Path to the folder with varying length RNAs in')

    parser_all = subparsers.add_parser('all', help='Run all of the commands one after the other')
    parser_all.add_argument('-a', '--adapter', help='Sequence of the adapter to remove from the 3\' end')
    parser_all.add_argument('-g', '--front', help='Sequence of the adapter to remove from the 5\' end')
    parser_all.add_argument('-b', '--anywhere', help='Sequence of the adapters to remove from both ends')
    parser_all.add_argument('-c', '--cutoff', help='Quality cutoff to trin RNA sequences at', default=20)
    parser_all.add_argument('-l', '--min-length', help='Minimum length to bin', type=int)
    parser_all.add_argument('-x', '--max-length', help='Maximum length to bin', type=int)
    parser_all.add_argument('-r', '--refseq', help='References for use with unitas', nargs='*')
    parser_all.add_argument('-u', '--subsiquent-refseq', help='Refseqs to use in a second run', nargs='*')
    parser_all.add_argument('-s', '--species', help='Species to set in unitas arguments', default='x')
    parser_all.add_argument('small_rna', help='Path to FASTQ containing the small RNA')
    parser_all.add_argument('genome', help='Genome to align against')

    args = parser.parse_args()

    load_config(args.path_to_config)

    def get_command_args(name):
        arguments = vars(args)

        if name in arguments and arguments[name] is not None:
            return arguments[name]
        else:
            try:
                return get_config_key('command', name)
            except KeyError:
                return None

    print(get_config_key('general', 'output_directory'))
    mkdir_if_not_exists(get_config_key('general', 'output_directory'))

    if args.command == 'process':
        process_command(
            get_command_args('small_rna'), 
            get_command_args('adapter'),
            get_command_args('front'),
            get_command_args('anywhere'),
            get_command_args('cutoff'),
            get_command_args('quiet')
        )

    if args.command == 'sort':
        sort_command(
            get_command_args('genome'),
            get_command_args('small_rna'),
            get_command_args('min_length'),
            get_command_args('max_length'),
            get_command_args('quiet')
        )

    if args.command == 'unitas':
        unitas_command(
            get_command_args('path_to_rnas'),
            get_command_args('species'),
            get_command_args('refseq'),
            get_command_args('subsiquent_refseq'),
            get_command_args('quiet')
        )

    if args.command == 'all':
        process_command(
            get_command_args('small_rna'), 
            get_command_args('adapter'),
            get_command_args('front'),
            get_command_args('anywhere'),
            get_command_args('cutoff'),
            get_command_args('quiet')
        )

        sort_command(
            get_command_args('genome'),
            os.path.join(get_config_key('general', 'output_directory'), 'cut_sequences.fastq'),
            get_command_args('min_length'),
            get_command_args('max_length'),
            get_command_args('quiet')
        )

        unitas_command(
            os.path.join(get_config_key('general', 'output_directory'), 'binned_rna'),
            get_command_args('species'),
            get_command_args('refseq'),
            get_command_args('subsiquent_refseq'),
            get_command_args('quiet')
        )
