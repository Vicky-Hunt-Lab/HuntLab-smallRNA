import os.path

from os import mkdir
from posixpath import commonpath
from subprocess import run

from trim import run_trim
from fastqc import run_fastqc, cut_rna_below_cutoff
from genome_align import align_to_genome, bin_rna_size

from config import get_config_key

from argparse import ArgumentParser

def process_command(small_rna, adapter, front, anywhere, cutoff, quiet):
    print('==> Starting command Process')
    small_rna_path = small_rna

    if adapter is not None or front is not None or anywhere is not None:
        run_trim(small_rna_path, adapter, front, anywhere)
        small_rna_path = os.path.join(get_config_key('general', 'output_directory'), 'trimmed_rna.fastq')
    else:
        print('==> No adapter sequence provided, skipping trim step')

    run_fastqc(small_rna_path, quiet=quiet)
    cut_rna_below_cutoff(small_rna_path, cutoff)

    print('==> Completed command Process')

def process_sort(genome, small_rna, quiet):
    print('==> Starting command Sort')
    new_fastq = align_to_genome(genome, small_rna, quiet=quiet)
    bin_rna_size(new_fastq)

    print('==> Completed command Sort')

if __name__ == '__main__':
    parser = ArgumentParser(description='stuff')
    parser.add_argument('-q', '--quiet', help='Supress output from intermediate commands', action='store_true')

    subparsers = parser.add_subparsers(dest='command')

    parser_process = subparsers.add_parser('process', help='Preprocessing for the RNA')
    parser_process.add_argument('-a', '--adapter', help='Sequence of the adapter to remove from the 3\' end')
    parser_process.add_argument('-g', '--front', help='Sequence of the adapter to remove from the 5\' end')
    parser_process.add_argument('-b', '--anywhere', help='Sequence of the adapters to remove from both ends')
    parser_process.add_argument('-c', '--cutoff', help='Quality cutoff to trin RNA sequences at', default=20, type=int)
    parser_process.add_argument('small_rna', help='Path to FASTQ containing the small RNA')

    parser_sort = subparsers.add_parser('sort', help='Find RNAs that align to a genome and sort them by length')
    parser_sort.add_argument('small_rna', help='Path to FASTQ containing the small RNA')
    parser_sort.add_argument('genome', help='Genome to align against')

    parser_classify = subparsers.add_parser('classify', help='Classify RNAs and run differntial expression')
    
    parser_multi = subparsers.add_parser('multi', help='Run multiple of the previous steps piping the output of one into the other')
    
    parser_all = subparsers.add_parser('all', help='Run all of the commands one after the other')
    parser_all.add_argument('-a', '--adapter', help='Sequence of the adapter to remove from the 3\' end')
    parser_all.add_argument('-g', '--front', help='Sequence of the adapter to remove from the 5\' end')
    parser_all.add_argument('-b', '--anywhere', help='Sequence of the adapters to remove from both ends')
    parser_all.add_argument('-c', '--cutoff', help='Quality cutoff to trin RNA sequences at', default=20)
    parser_all.add_argument('small_rna', help='Path to FASTQ containing the small RNA')
    parser_all.add_argument('genome', help='Genome to align against')

    args = parser.parse_args()

    def get_command_args(name):
        arguments = vars(args)

        if name in arguments and arguments[name] is not None:
            return arguments[name]
        else:
            try:
                return get_config_key('command', name)
            except KeyError:
                return None

    mkdir(get_config_key('general', 'output_directory'))

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
        process_sort(
            get_command_args('genome'),
            get_command_args('small_rna'),
            get_command_args('quiet')
        )

    if args.command == 'classify':
        pass

    if args.command == 'all':
        process_command(
            get_command_args('small_rna'), 
            get_command_args('adapter'),
            get_command_args('front'),
            get_command_args('anywhere'),
            get_command_args('cutoff'),
            get_command_args('quiet')
        )

        process_sort(
            get_command_args('genome'),
            os.path.join(get_config_key('general', 'output_directory'), 'cut_sequences.fastq'),
            get_command_args('quiet')
        )