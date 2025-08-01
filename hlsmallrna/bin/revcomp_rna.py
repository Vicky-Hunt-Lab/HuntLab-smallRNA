#!/usr/bin/env python3
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

from argparse import ArgumentParser

from Bio import SeqIO

def reverse_complement(seq):
    revcomp = seq.reverse_complement()
    revcomp.id = seq.id
    revcomp.description = seq.description

    return revcomp

def main():
    parser = ArgumentParser(description='Simple script to reverse complement a FASTA file')

    parser.add_argument('-i', '--input', help='FASTA file to reverse complement')
    parser.add_argument('-o', '--output', help='FASTA file to write to')

    args = parser.parse_args()

    if args.input is None or args.output is None:
        raise Exception('Both -i and -o need to be set to run this')
        
    seqs = SeqIO.parse(args.input, 'fasta')
    seqs = map(reverse_complement, seqs)

    SeqIO.write(seqs, args.output, 'fasta')
