# Copyright 2022 Vicky Hunt Lab Members
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#!/usr/bin/env python3

import os.path
import glob
import shutil

from math import inf
from argparse import ArgumentParser

from Bio import SeqIO, SeqRecord
from gffutils.iterators import DataIterator

from .trim import run_trim
from .fastqc import run_fastqc, cut_rna_below_cutoff
from .genome_align import align_to_genome, bin_rna_size, graph_length
from .create_noncoding import extract_noncoding
from .unitas import run_unitas_annotation, merge_summary, graph_unitas_classification_type
from .targetid import revcomp_input_file, find_targets, build_summary_files
from .ss_overlap import samestrand_overlap

from .config import get_config_key, mkdir_if_not_exists, load_config, do_log

def validate_file(path_to_file, file_type):
    '''
    Examine a file and determine if it is of the requested file type
    '''
    if file_type == 'directory':
        return os.path.isdir(path_to_file)
    elif file_type == 'fasta':
        try:
            item = next(SeqIO.parse(path_to_file, 'fasta'))
            return type(item) == SeqRecord.SeqRecord
        except:
            return False
    elif file_type == 'fastq':
        try:
            item = next(SeqIO.parse(path_to_file, 'fastq'))
            return type(item) == SeqRecord.SeqRecord
        except:
            return False
    elif file_type == 'gff':
        iter = DataIterator(path_to_file)
        return len(iter.peek) > 0


def process_command(small_rna, adapter, front, anywhere, cutoff, quiet):
    '''
    Code to run when the user chooses the process command
    '''

    if not validate_file(small_rna, 'fastq'):
        print(f'Error: expected a small RNA FASTQ with at least one sequence, got {small_rna}')
        return False

    do_log(quiet, '==> Starting command Process')
    small_rna_path = small_rna

    if adapter is not None or front is not None or anywhere is not None:
        run_trim(small_rna_path, adapter, front, anywhere)
        small_rna_path = os.path.join(get_config_key('general', 'output_directory'), 'trimmed_rna.fastq')
    else:
        do_log(quiet, '==> No adapter sequence provided, skipping trim step')

    if cutoff > 0:
        run_fastqc(small_rna_path, quiet=quiet)
        cut_rna_below_cutoff(small_rna_path, cutoff, quiet=quiet)
    else:
        # create the cut_sequences.fastq, even if fastqc wasn't run
        shutil.copy2(small_rna_path, os.path.join(get_config_key('general', 'output_directory'), 'cut_sequences.fastq'))

    do_log(quiet, '==> Completed command Process')

def sort_command(genome, small_rna, cds, min_length, max_length, quiet):
    '''
    Code to run when the user chooses the sort command
    '''

    if not validate_file(genome, 'fasta'):
        print(f'Error: expected a genome in FASTA format, got {genome}')
        return False

    if not validate_file(small_rna, 'fastq'):
        print(f'Error: expected a small RNA FASTQ with at least one sequence, got {small_rna}')
        return False

    do_log(quiet, '==> Starting command Sort')
    new_fastq = align_to_genome(genome, small_rna, cds, quiet=quiet)
    table_file = bin_rna_size(new_fastq, min_length, max_length, quiet=quiet)

    graph_length(table_file)

    do_log(quiet, '==> Completed command Sort')

def extractnc_command(genome, gff, quiet):
    '''
    Code to run when the user chooses to extract the noncoding mRNA reigon
    '''

    if not validate_file(genome, 'fasta'):
        print(f'Error: expected a genome in FASTA format, got {genome}')
        return False

    if not validate_file(gff, 'gff'):
        print(f'Error: expected a GFF file with at least one feature, got {gff}')
        return False

    do_log(quiet, '==> Starting command extractNC')

    extract_noncoding(
        genome, gff, 
        quiet=quiet,
        output=os.path.join(get_config_key('general', 'output_directory'), 'noncoding.fasta')
    )

    do_log(quiet, '==> Completed command extractNC')

def unitas_command(small_rna_path, species_name, ref_seqs, quiet):
    '''
    Code to run when the user chooses the unitas command
    '''

    if not validate_file(small_rna_path, 'directory'):
        print(f'Error: expected a directory containing small RNA FASTQ of varing lengths, got {small_rna_path}')
        return False

    if species_name == 'x' and ref_seqs is None:
        print(f'Error: expected at least one of, a non x species name (-s) or at least one reference file (-r)')
        return False

    do_log(quiet, '==> Starting command Unitas')
    UNITAS_OUTPUT = os.path.join(get_config_key('general', 'output_directory'), 'unitas')

    mkdir_if_not_exists(UNITAS_OUTPUT)

    for small_rna in glob.glob(os.path.join(small_rna_path, '*.fastq')):
        run_unitas_annotation(small_rna, species_name, ref_seqs, quiet=quiet, unitas_output=UNITAS_OUTPUT)

    table_path = merge_summary()
    graph_unitas_classification_type(table_path)
    do_log(quiet, '==> Completed command Unitas')

def targetid_command(small_rna, targets, min_seq_length, mismatches_allowed, quiet):
    '''
    Code to run when the user chooses the targetid command
    '''

    if not validate_file(small_rna, 'fastq'):
        print(f'Error: expected a small RNA FASTQ with at least one sequence, got {small_rna}')
        return False

    if targets is None:
        print(f'Error: expected at least one target file (-t)')
        return False

    do_log(quiet, '==> Starting TargetID command')

    if targets is None:
        print('Error: You need to supply at least one target file with -t')

    revcomp_file = revcomp_input_file(small_rna, quiet=quiet)
    sam_files = find_targets(revcomp_file, targets, min_seq_length=min_seq_length, mismatches_allowed=mismatches_allowed, quiet=quiet)
    build_summary_files(sam_files, quiet=quiet)

    do_log(quiet, '==> Ending TargetID command')

def main():
    parser = ArgumentParser(description='Pipeline to process small RNAs')
    parser.add_argument('-q', '--quiet', help='Supress output from intermediate commands', action='count', default=0)
    parser.add_argument('-C', '--path-to-config', help='Path to the TOML format config file to use', default='config.toml')

    subparsers = parser.add_subparsers(dest='command')

    parser_process = subparsers.add_parser('process', help='Preprocessing for the RNA')
    parser_process.add_argument('-a', '--adapter', help='Sequence of the adapter to remove from the 3\' end')
    parser_process.add_argument('-g', '--front', help='Sequence of the adapter to remove from the 5\' end')
    parser_process.add_argument('-b', '--anywhere', help='Sequence of the adapters to remove from both ends')
    parser_process.add_argument('-c', '--cutoff', help='Quality cutoff to trin RNA sequences at', default=20, type=int)
    parser_process.add_argument('small_rna', help='Path to FASTQ containing the small RNA')

    parser_sort = subparsers.add_parser('sort', help='Find RNAs that align to a genome and sort them by length')
    parser_sort.add_argument('-d', '--cds', help='Optional CDS region, also align this to the CDS reigon as well as the genome')
    parser_sort.add_argument('-l', '--min-length', help='Minimum length to bin', type=int, default=-inf)
    parser_sort.add_argument('-x', '--max-length', help='Maximum length to bin', type=int, default=inf)
    parser_sort.add_argument('small_rna', help='Path to FASTQ containing the small RNA')
    parser_sort.add_argument('genome', help='Genome to align against')

    parser_extractnc = subparsers.add_parser('extractnc', help='Extarct the noncoding reigon from a fasta with a GFF file')
    parser_extractnc.add_argument('genome', help='FASTA containing the genome to extract from')
    parser_extractnc.add_argument('gff_file', help='GFF file containing annotations of CDS and mRNA reigons')

    parser_unitas = subparsers.add_parser('unitas', help='Run unitas on split files and merge results')
    parser_unitas.add_argument('-r', '--refseq', help='References for use with unitas', nargs='*', default=None)
    parser_unitas.add_argument('-s', '--species', help='Species to set in unitas arguments', default='x')
    parser_unitas.add_argument('path_to_rnas', help='Path to the folder with varying length RNAs in')

    parser_targetid = subparsers.add_parser('targetid', help='Align small RNA to a number of genome features to find out what is targeted')
    parser_targetid.add_argument('-m', '--min-seq-length', help='Minimum sequence length to properly align', default=5)
    parser_targetid.add_argument('-t', '--target-files', help='Files containing genome features that could be targeted', nargs='+')
    parser_targetid.add_argument('--num-mismatches', help='Number of mismatches to allow in the alignment, defaults to 0', default=0, type=int)
    parser_targetid.add_argument('small_rna', help='Path to the FASTQ containing the small RNA to find targets of')

    parser_all = subparsers.add_parser('all', help='Run process, sort and unitas one after the other')
    parser_all.add_argument('-a', '--adapter', help='Sequence of the adapter to remove from the 3\' end')
    parser_all.add_argument('-d', '--cds', help='Optional CDS region, also align this to the CDS reigon as well as the genome')
    parser_all.add_argument('-g', '--front', help='Sequence of the adapter to remove from the 5\' end')
    parser_all.add_argument('-b', '--anywhere', help='Sequence of the adapters to remove from both ends')
    parser_all.add_argument('-c', '--cutoff', help='Quality cutoff to trin RNA sequences at', default=20, type=int)
    parser_all.add_argument('-l', '--min-length', help='Minimum length to bin', type=int, default=-inf)
    parser_all.add_argument('-x', '--max-length', help='Maximum length to bin', type=int, default=inf)
    parser_all.add_argument('-r', '--refseq', help='References for use with unitas', nargs='*', default=None)
    parser_all.add_argument('-s', '--species', help='Species to set in unitas arguments', default='x')
    parser_all.add_argument('small_rna', help='Path to FASTQ containing the small RNA')
    parser_all.add_argument('genome', help='Genome to align against')

    args = parser.parse_args()

    do_log(args.quiet, f"Output directory set to {os.path.abspath(get_config_key('general', 'output_directory'))}")
    if os.path.exists(get_config_key('general', 'output_directory')):
        do_log(args.quiet, 'WARNING: output directory already exists, this run may use previously generated data, if you didn\'t mean this, try again with a different output directory')

    load_config(args.path_to_config, quiet=args.quiet)

    def get_command_args(name):
        arguments = vars(args)

        if name in arguments and arguments[name] is not None:
            return arguments[name]
        else:
            try:
                return get_config_key('command', name)
            except KeyError:
                return None

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
            get_command_args('cds'),
            get_command_args('min_length'),
            get_command_args('max_length'),
            get_command_args('quiet')
        )

    if args.command == 'extractnc':
        extractnc_command(
            get_command_args('genome'),
            get_command_args('gff_file'),
            get_command_args('quiet')
        )

    if args.command == 'unitas':
        unitas_command(
            get_command_args('path_to_rnas'),
            get_command_args('species'),
            get_command_args('refseq'),
            get_command_args('quiet')
        )

    if args.command == 'targetid':
        targetid_command(
            get_command_args('small_rna'),
            get_command_args('target_files'),
            get_command_args('min_seq_length'),
            get_command_args('num_mismatches'),
            get_command_args('quiet')
        )

    if args.command == 'all':
        out_code = process_command(
            get_command_args('small_rna'), 
            get_command_args('adapter'),
            get_command_args('front'),
            get_command_args('anywhere'),
            get_command_args('cutoff'),
            get_command_args('quiet')
        )

        if out_code is not None:
            return

        out_code = sort_command(
            get_command_args('genome'),
            os.path.join(get_config_key('general', 'output_directory'), 'cut_sequences.fastq'),
            get_command_args('cds'),
            get_command_args('min_length'),
            get_command_args('max_length'),
            get_command_args('quiet')
        )

        if out_code is not None:
            return

        unitas_command(
            os.path.join(get_config_key('general', 'output_directory'), 'binned_rna'),
            get_command_args('species'),
            get_command_args('refseq'),
            get_command_args('quiet')
        )

def main_ssoverlap():
    '''
    Main function form the overlap_ss script
    '''
    
    parser = ArgumentParser(description='Looks for overlaps of two classes of RNA on the same strand of DNA')

    parser.add_argument('genome', help='FASTA file containing the genome of the organism')
    parser.add_argument('rna_file_1', help='File containing the RNA to use as a base')
    parser.add_argument('rna_file_2', help='File containing the RNA to look for overlaps with')

    parser.add_argument('-q', '--quiet', help='Supress output from intermediate commands', action='count', default=0)    

    args = parser.parse_args()

    OUTPUT_DIR = os.path.join(get_config_key('general', 'output_directory'), 'samestrand_overlap')

    mkdir_if_not_exists(OUTPUT_DIR)

    do_log(args.quiet, '==> Converting RNA to FASTA format')

    global longest_rna
    longest_rna = 0
    for rna_file in [args.rna_file_1, args.rna_file_2]:
        if (
            rna_file.endswith('.fa') 
            or rna_file.endswith('.fasta') 
            or rna_file.endswith('.fna') 
            or rna_file.endswith('.ffn') 
            or rna_file.endswith('.faa') 
            or rna_file.endswith('.frn')
        ):
            filetype = 'fasta'
        elif rna_file.endswith('.fq') or rna_file.endswith('.fastq'):
            filetype = 'fastq'

        def convert_to_fasta(seq):
            global longest_rna
            if len(seq) > longest_rna:
                longest_rna = len(seq)

            return seq

        seqs = SeqIO.parse(rna_file, filetype)
        seqs = map(convert_to_fasta, seqs)
        SeqIO.write(seqs, os.path.join(OUTPUT_DIR, os.path.basename(rna_file) + '.fasta'), 'fasta')

    samestrand_overlap(
        args.genome, 
        os.path.join(OUTPUT_DIR, os.path.basename(args.rna_file_1) + '.fasta'), 
        os.path.join(OUTPUT_DIR, os.path.basename(args.rna_file_2) + '.fasta'),
        longest_rna,
        quiet=args.quiet
    )


if __name__ == '__main__':
    main()