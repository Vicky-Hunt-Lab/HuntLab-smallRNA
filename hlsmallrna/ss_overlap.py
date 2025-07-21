# Copyright 2022 - 2025 Vicky Hunt Lab Members
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
import os

from subprocess import run
from argparse import ArgumentParser

from Bio import SeqIO

def samestrand_overlap(genome_file, rna_file_1, rna_file_2, longest_rna, verbose=False):
    '''
    Run the samestand overlap script
    '''
    print('==> Running same strand overlap script')

    CWD = os.getcwd()
    OUTPUT_DIR = os.path.join(CWD, 'samestrand_overlap_output')

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    os.chdir(OUTPUT_DIR)

    command = [
        'overlap_ss.sh',
        os.path.join(CWD, genome_file),
        os.path.join(CWD, rna_file_1),
        os.path.join(CWD, rna_file_2),
        str(longest_rna)
    ]

    run(command, capture_output=not verbose)

    os.chdir(CWD)

    print('==> Finished same strand overlap script')

def main_ssoverlap():
    '''
    Main function form the overlap_ss script
    '''
    
    parser = ArgumentParser(description='Looks for overlaps of two classes of RNA on the same strand of DNA')

    parser.add_argument('genome', help='FASTA file containing the genome of the organism')
    parser.add_argument('rna_file_1', help='File containing the RNA to use as a base')
    parser.add_argument('rna_file_2', help='File containing the RNA to look for overlaps with')

    parser.add_argument('-v', '--verbose', action='store_true', help='Print all output for this command')    
    parser.add_argument('-o', '--output', help='Directory to output to', default='samestrand_overlap')

    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    print('==> Converting RNA to FASTA format')

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
        SeqIO.write(seqs, os.path.join(args.output, os.path.basename(rna_file) + '.fasta'), 'fasta')

    samestrand_overlap(
        args.genome, 
        os.path.join(args.output, os.path.basename(args.rna_file_1) + '.fasta'), 
        os.path.join(args.output, os.path.basename(args.rna_file_2) + '.fasta'),
        longest_rna,
        verbose=args.verbose
    )