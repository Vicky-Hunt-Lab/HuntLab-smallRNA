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

import csv
from os import chroot, path
import re

from argparse import ArgumentParser
from collections import defaultdict

import pysam

from Bio import SeqIO
from gffutils.iterators import DataIterator

def process_chrom_biopython_merge_all(path_to_fasta, filetype, output, band_width=20000):
    '''
    Use Biopython to produce the genome coordinates, version that merges all scaffolds into one 
    '''
    sequences = SeqIO.parse(path_to_fasta, filetype)

    scaffolds = {}

    result_file = []
    scaffold_info = []
    for seq in sequences:
        scaffolds[seq.id] = {'name': seq.id, 'length': len(seq)}

    current_start = 1
    for scaffold in sorted(scaffolds, key=lambda x: x['length']):
        result_file.append([
                'GENOME',
                current_start, current_start + scaffolds[scaffold]['length'],
                scaffold
            ])

        current_start = current_start + scaffolds[scaffold]['length'] + 1
        scaffold_info.append(['GENOME', current_start, current_start + band_width, scaffold, 'red'])

    with open(output, 'w') as f:
        writer = csv.writer(f, delimiter='\t')

        writer.writerow(['Chrom', 'Start', 'End', 'Name'])
        writer.writerows(result_file)
        
    with open('scaffold_info_' + output, 'w') as f:
        writer = csv.writer(f, delimiter='\t')

        writer.writerow(['Chrom', 'Start', 'End', 'Name', 'Colors'])
        writer.writerows(scaffold_info)

def process_chrom_biopython(path_to_fasta, filetype, output, band_width=20000, scaffold_aware=False):
    '''
    Use Biopython to produce the genome coordinates, scaffold aware tries to merge scaffolds into continuious reigons
    '''
    sequences = SeqIO.parse(path_to_fasta, filetype)

    scaffolds = defaultdict(lambda: {})

    result_file = []
    scaffold_info = []
    for seq in sequences:
        matches = re.findall(r'(.+)_scaffold(.+)', seq.id)
        if not scaffold_aware or len(matches) == 0:
            result_file.append([seq.id, 1, len(seq), seq.id])

            scaffold_info.append([seq.id, 1, len(seq), seq.id, 'grey'])
        else:
            scaffolds[matches[0][0]][int(matches[0][1])] = {'name': seq.id, 'length': len(seq)}

    for extra in scaffolds:
        current_start = 1

        for key in sorted(scaffolds[extra].keys()):
            result_file.append([
                extra,
                current_start, current_start + scaffolds[extra][key]['length'],
                extra + f'_scaffold{key}'
            ])

            current_start = current_start + scaffolds[extra][key]['length'] + 1
            scaffold_info.append([extra, current_start, current_start + band_width, extra + f'_scaffold{key}', 'red'])
            
    with open(output, 'w') as f:
        writer = csv.writer(f, delimiter='\t')

        writer.writerow(['Chrom', 'Start', 'End', 'Name'])
        writer.writerows(result_file)

    if len(scaffold_info) > 0:
        with open('scaffold_info_' + output, 'w') as f:
            writer = csv.writer(f, delimiter='\t')

            writer.writerow(['Chrom', 'Start', 'End', 'Name', 'Colors'])
            writer.writerows(scaffold_info)
        

def process_feature_pysam(path_to_sam, filetype, output, scaffold_index=None):
    '''
    Use pysam to produce coordinate files from SAM and BAM files
    '''
    read_flag = 'r'
    
    if filetype == 'bam':
        read_flag += 'b'

    data = pysam.AlignmentFile(path_to_sam, read_flag)
    fileheader = data.header
    
    sense_targets = []
    antisene_targets = []

    for alignment in data:
        if not alignment.is_unmapped:
            if scaffold_index is None:
                if alignment.is_reverse:
                    antisene_targets.append([
                        fileheader.get_reference_name(alignment.reference_id), alignment.query_name, 
                        alignment.reference_start + 1, alignment.reference_start + alignment.query_alignment_end + 1
                    ])
                else:
                    sense_targets.append([
                        fileheader.get_reference_name(alignment.reference_id), alignment.query_name, 
                        alignment.reference_start + 1, alignment.reference_start + alignment.query_alignment_end + 1
                    ])
            else:
                name = fileheader.get_reference_name(alignment.reference_id)
                start = alignment.reference_start + 1
                end = alignment.reference_start + alignment.query_alignment_end + 1

                if name in scaffold_index:
                    start += scaffold_index[name][1]
                    end += scaffold_index[name][1]

                    name = scaffold_index[name][0]

                if alignment.is_reverse:
                    antisene_targets.append([
                        name, alignment.query_name, 
                        start, end
                    ])
                else:
                    sense_targets.append([
                        name, alignment.query_name, 
                        start, end
                    ])

    with open('sense_' + output, 'w') as f:
        writer = csv.writer(f, delimiter='\t')

        writer.writerow(['Chrom', 'Name', 'Start', 'End'])
        writer.writerows(sense_targets)

    with open('antisense_' + output, 'w') as f:
        writer = csv.writer(f, delimiter='\t')

        writer.writerow(['Chrom', 'Name', 'Start', 'End'])
        writer.writerows(antisene_targets)

def process_feature_gffutils(path_to_gff, output, feature, scaffold_index=None):
    '''
    Use gffutils to produce coordinate files from GFF files
    '''
    data = DataIterator(path_to_gff)
    
    sense_targets = []
    antisene_targets = []

    for i, line in enumerate(data):
        if feature is None or line[2] == feature:
            if scaffold_index is None:
                if line[6] == '+':
                    sense_targets.append([line[0], line[2], line[3], line[4]])
                elif line[6] == '-':
                    antisene_targets.append([line[0], line[2], line[3], line[4]])
            else:
                name = line[0]
                start = int(line[3])
                end = int(line[4])

                if name in scaffold_index:
                    start += scaffold_index[name][1]
                    end += scaffold_index[name][1]

                    name = scaffold_index[name][0]

                if line[6] == '-':
                    antisene_targets.append([
                        name, line[2],
                        start, end
                    ])
                elif line[6] == '+':
                    sense_targets.append([
                        name, line[2],
                        start, end
                    ])

    with open('sense_' + output, 'w') as f:
        writer = csv.writer(f, delimiter='\t')

        writer.writerow(['Chrom', 'Name', 'Start', 'End'])
        writer.writerows(sense_targets)

    with open('antisense_' + output, 'w') as f:
        writer = csv.writer(f, delimiter='\t')

        writer.writerow(['Chrom', 'Name', 'Start', 'End'])
        writer.writerows(antisene_targets)

def process_fa_out(path_to_fa_out, output, feature, scaffold_index=None):
    '''
    Custom function to parse RepeatMasker .fa.out files
    '''

    sense_targets = []
    antisense_targets = []

    with open(path_to_fa_out) as f:
        for i, line in enumerate(f):
            if i < 3:
                continue

            line = line.split()
            name = line[10]

            if scaffold_index is None:
                if line[8].lower() == 'c':
                    antisense_targets.append([
                        line[4], line[10],
                        line[5], line[6]
                    ])
                else:
                    sense_targets.append([
                        line[4], line[10],
                        line[5], line[6]
                    ])
            else:
                chrom = line[4]
                start = int(line[5])
                end = int(line[6])

                if chrom in scaffold_index:
                    start += scaffold_index[chrom][1]
                    end += scaffold_index[chrom][1]

                    chrom = scaffold_index[chrom][0]

                    if line[8].lower() == 'c':
                        antisense_targets.append([
                            chrom, line[10],
                            start, end
                        ])
                    else:
                        sense_targets.append([
                            chrom, line[10],
                            start, end
                        ])

    with open('sense_' + output, 'w') as f:
        writer = csv.writer(f, delimiter='\t')

        writer.writerow(['Chrom', 'Name', 'Start', 'End'])
        writer.writerows(sense_targets)

    with open('antisense_' + output, 'w') as f:
        writer = csv.writer(f, delimiter='\t')

        writer.writerow(['Chrom', 'Name', 'Start', 'End'])
        writer.writerows(antisense_targets)

def main():
    parser = ArgumentParser(description='Convert common bioinformatics file formats into coordinate file to plot')

    parser.add_argument('-a', '--fasta', action='store_true', help='Treat input as a FASTA file')
    parser.add_argument('-q', '--fastq', action='store_true', help='Treat input as a FASTQ file')
    parser.add_argument('-s', '--sam', action='store_true', help='Treat input as a SAM file')
    parser.add_argument('-b', '--bam', action='store_true', help='Treat input as a BAM file')
    parser.add_argument('-g', '--gff', action='store_true', help='Treat input as a GFF file')
    parser.add_argument('-r', '--rm-fa-out', action='store_true', help='Treat input as a RepeatMasker .fa.out file')

    parser.add_argument('-c', '--scaffold-aware', action='store_true', help='Merge scaffolds into one chromosome')
    parser.add_argument('-m', '--merge-all', action='store_true', help='When building the genome file merge all scaffolds')

    parser.add_argument('--feature', help='Select a gff feature to use')
    parser.add_argument('--genome-coords', help='File containing the coordinates of the genome')
    parser.add_argument('--band-width', help='Size of the bands showing change in scaffolds', type=int, default=20000)

    parser.add_argument('-o', '--output', help='File to output to', default='result.tsv')
    parser.add_argument('input_file', help='File to convert, attempts to autodetect type')

    args = parser.parse_args()

    filetype = None
    if args.fasta:
        filetype = 'fasta'
    elif args.fastq:
        filetype = 'fastq'
    elif args.sam:
        filetype = 'sam'
    elif args.bam:
        filetype = 'bam'
    elif args.gff:
        filetype = 'gff'
    elif args.rm_fa_out:
        filetype = 'fa_out'

    if filetype is None:
        if (
            args.input_file.endswith('.fa') 
            or args.input_file.endswith('.fasta') 
            or args.input_file.endswith('.fna') 
            or args.input_file.endswith('.ffn') 
            or args.input_file.endswith('.faa') 
            or args.input_file.endswith('.frn')
        ):
            filetype = 'fasta'
        elif args.input_file.endswith('.fq') or args.input_file.endswith('.fastq'):
            filetype = 'fastq'
        elif args.input_file.endswith('.sam'):
            filetype = 'sam'
        elif args.input_file.endswith('.bam'):
            filetype = 'bam'
        elif args.input_file.endswith('.gff') or args.input_file.endswith('.gff3'):
            filetype = 'gff'
        elif args.input_file.endswith('.fa.out'):
            filetype = 'fa_out'

    if filetype is None:
        raise Exception('Filetype could not be determined, try explicitly setting an option')

    print(f'File type detected as : {filetype}')

    if filetype in ['fasta', 'fastq']:
        if not args.merge_all:
            process_chrom_biopython(args.input_file, filetype, args.output, band_width=args.band_width, scaffold_aware=args.scaffold_aware)
        else:
            process_chrom_biopython_merge_all(args.input_file, filetype, args.output, band_width=args.band_width)
    else:
        scaffold_index = None

        if args.scaffold_aware and args.genome_coords is None:
            raise Exception('For non-fasta file types, you need to set --genome-coords if you set --scaffold-aware')
        elif args.scaffold_aware:
            with open(args.genome_coords) as f:
                reader = csv.reader(f, delimiter='\t')
                scaffold_index = {}

                next(reader)
                for line in reader:
                    scaffold_index[line[-1]] = (line[0], int(line[1]))

        if filetype in ['sam', 'bam']:
            process_feature_pysam(args.input_file, filetype, args.output, scaffold_index=scaffold_index)
        elif filetype == 'gff':
            process_feature_gffutils(args.input_file, args.output, args.feature, scaffold_index=scaffold_index)
        elif filetype == 'fa_out':
            process_fa_out(args.input_file, args.output, args.feature, scaffold_index=scaffold_index)
