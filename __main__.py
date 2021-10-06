import os.path

from os import mkdir
from posixpath import commonpath
from subprocess import run

from trim import run_trim
from fastqc import run_fastqc, cut_rna_below_cutoff

from config import get_config_key

from argparse import ArgumentParser

if __name__ == '__main__':
    parser = ArgumentParser(description='stuff')

    subparsers = parser.add_subparsers(dest='command')

    parser_process = subparsers.add_parser('process', help='Preprocessing for the RNA')
    parser_process.add_argument('-a', '--adapter', help='Sequence of the adapter to remove from the 3\' end')
    parser_process.add_argument('-g', '--front', help='Sequence of the adapter to remove from the 5\' end')
    parser_process.add_argument('-b', '--anywhere', help='Sequence of the adapters to remove from both ends')
    parser_process.add_argument('-c', '--cutoff', help='Quality cutoff to trin RNA sequences at', default=20)
    parser_process.add_argument('small_rna', help='Path to FASTQ containing the small RNA')

    parser_sort = subparsers.add_parser('sort', help='Find RNAs that align to a genome and sort them by length')

    parser_classify = subparsers.add_parser('classify', help='Classify RNAs and run differntial expression')
    
    parser_multi = subparsers.add_parser('multi', help='Run multiple of the previous steps piping the output of one into the other')
    
    parser_all = subparsers.add_parser('all', help='Run all of the commands one after the other')

    args = parser.parse_args()

    mkdir(get_config_key('general', 'output_directory'))

    if args.command == 'process':
        small_rna_path = args.small_rna

        if args.adapter is not None or args.front is not None or args.anywhere is not None:
            run_trim(small_rna_path, args.adapter, args.front, args.anywhere)
            small_rna_path = os.path.join(get_config_key('general', 'output_directory'), 'trimmed_rna.fastq')

        run_fastqc(small_rna_path)
        cut_rna_below_cutoff(small_rna_path, args.cutoff)


    