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

def prepend_label_to_id(label, seqs):
    '''
    Add the specified label to the start of the sequence ID
    '''
    for seq in seqs:
        seq.id = label + '|' + seq.id

        yield seq

def label_file_for_unitas(input, output, label):
    '''
    Prepend a specified label to all sequences in a FASTA file
    '''
    SeqIO.write(
        prepend_label_to_id(label, SeqIO.parse(input, 'fasta')),
        output,
        'fasta'
    )

def label_for_unitas_cli():
    parser = ArgumentParser(description='Program to automatically prepend unitas labels to all the sequences in a FASTA file')

    parser.add_argument('-o', '--output', help='File to output to, defaults to labelled.fasta', default='labeled.fasta')

    parser.add_argument('label', help='Label to prepend to the sequences in a fasta file')
    parser.add_argument('file_path', help='Path to the fasta file to label')

    args = parser.parse_args()

    label_file_for_unitas(args.file_path, args.output, args.label)