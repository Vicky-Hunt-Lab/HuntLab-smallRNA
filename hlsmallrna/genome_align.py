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
from shutil import copy2

import os
import os.path
import csv
import shutil

from subprocess import run
from collections import defaultdict

import numpy as np

from Bio import SeqIO, SeqRecord, Seq
from matplotlib import pyplot as plt

def into_seqrecord(seqs):
    '''
    Convert a set of Biopython sequences into SeqRecord objects
    '''
    for i, seq in enumerate(seqs):
        record = SeqRecord.SeqRecord(seq=Seq.Seq(seq), id=f'SmallRNA{i:07}', description='')
        record.letter_annotations["phred_quality"] = [40] * len(record)

        yield record

def merge_genome_cds_hits(results_seq, cds_seq):
    '''
    Merge FASTQ hits from the genome and CDS
    '''
    for seq in SeqIO.parse(results_seq, 'fastq'):
        yield seq.seq

    if os.path.exists(cds_seq):
        for seq in SeqIO.parse(cds_seq, 'fastq'):
            yield seq.seq

def align_to_genome(genome, small_rnas, cds, program_paths, output, verbose=False, keep_intermediate=False, threads=4, mismatches=0):
    '''
    Align the small RNAs to the genome and filter out any that are unsuccessful
    '''

    INTERMEDIATE_SAM = os.path.join(output, 'mapped_sequences.sam')
    INTERMEDIATE_BAM = os.path.join(output, 'mapped_sequences.bam')
    UNMAPPED_BAM = os.path.join(output, 'unmapped_sequences.bam')
    RESULT_FASTQ = os.path.join(output, 'genome_mapped_sequences.fastq')
    RESULT_UNMAPPED_FASTQ = os.path.join(output, 'genome_unmapped_sequences.fastq')
    INDEX_DIRECTORY = os.path.join(output, 'bowtie2_index')

    CDS_INTERMEDIATE_SAM = os.path.join(output, 'cds_sequences.sam')
    CDS_INTERMEDIATE_BAM = os.path.join(output, 'cds_sequences.bam')
    CDS_UNMAPPED_BAM = os.path.join(output, 'cds_unmapped_sequences.bam')
    CDS_UNMAPPED_FASTQ = os.path.join(output, 'cds_unmapped_sequences.fastq')
    CDS_RESULT_FASTQ = os.path.join(output, 'cds_sequences.fastq')

    FINAL_FASTQ = os.path.join(output, 'mapped_sequences.fastq')
    FINAL_UNMAPPED_FASTQ = os.path.join(output, 'unmapped_sequences.fastq')

    if not os.path.exists(INDEX_DIRECTORY):
        os.makedirs(INDEX_DIRECTORY)

    bowtie2_build_index = [
        program_paths['bowtie2-build'],
        '--threads', str(threads),
        genome,
        os.path.join(INDEX_DIRECTORY, 'genome_index')
    ]

    cds_build_index = [
        program_paths['bowtie2-build'],
        '--threads', str(threads),
        cds,
        os.path.join(INDEX_DIRECTORY, 'cds_index')
    ]

    if mismatches == 0:
        bowtie2_align_reads = [
                program_paths['bowtie2'],
                '--threads', str(threads),
                '-L', '18',
                '--no-1mm-upfront',
                '--score-min', 'L,0,0',
                '--end-to-end',
                '-M', '0',
                '-x', os.path.join(INDEX_DIRECTORY, 'genome_index'),
                '-U', small_rnas,
                '-S', INTERMEDIATE_SAM
            ]
    else:
        bowtie2_align_reads = [
            program_paths['bowtie2'],
            '--threads', str(threads),
            '-L', '18',
            '--score-min', 'L,' + str(-mismatches) + ',0',
            '--end-to-end',
            '--mp', '1,1',
            '--ignore-quals',
            '--rdg', '100,100',
            '--rfg', '100,100',
            '-x', os.path.join(INDEX_DIRECTORY, 'genome_index'),
            '-U', small_rnas,
            '-S', INTERMEDIATE_SAM
        ]

    if mismatches == 0:
        cds_align_reads = [
                program_paths['bowtie2'],
                '--threads', str(threads),
                '-L', '18',
                '--no-1mm-upfront',
                '--score-min', 'L,0,0',
                '--end-to-end',
                '-M', '0',
                '-x', os.path.join(INDEX_DIRECTORY, 'cds_index'),
                '-U', RESULT_UNMAPPED_FASTQ,
                '-S', CDS_INTERMEDIATE_SAM
            ]
    else:
        cds_align_reads = [
            program_paths['bowtie2'],
            '--threads', str(threads),
            '-L', '18',
            '--score-min', 'L,' + str(-mismatches) + ',0',
            '--end-to-end',
            '--mp', '1,1',
            '--ignore-quals',
            '--rdg', '100,100',
            '--rfg', '100,100',
            '-x', os.path.join(INDEX_DIRECTORY, 'cds_index'),
            '-U', RESULT_UNMAPPED_FASTQ,
            '-S', CDS_INTERMEDIATE_SAM
        ]

    samtools_view = [
        program_paths['samtools'],
        'view',
        '-h',
        '-b',
        '-F', '4',
        '-F', '256',
        '-o', INTERMEDIATE_BAM,
        INTERMEDIATE_SAM,
    ]

    bedtools_bamtofastq_command = [
        program_paths['samtools'],
        'fastq',
        INTERMEDIATE_BAM,
        '-o', RESULT_FASTQ,
        '-0', RESULT_FASTQ
    ]

    cds_samtools_view = [
        program_paths['samtools'],
        'view',
        '-h',
        '-b',
        '-F', '4',
        '-F', '256',
        '-o', CDS_INTERMEDIATE_BAM,
        CDS_INTERMEDIATE_SAM
    ]

    cds_bedtools_bamtofastq_command = [
        program_paths['samtools'],
        'fastq',
        CDS_INTERMEDIATE_BAM,
        '-o', CDS_RESULT_FASTQ,
        '-0', CDS_RESULT_FASTQ
    ]
    
    unmapped_samtools_view = [
        program_paths['samtools'],
        'view',
        '-h',
        '-b',
        '-f', '4',
        '-o', UNMAPPED_BAM,
        INTERMEDIATE_SAM
    ]

    unmapped_bedtools_bamtofastq_command = [
        program_paths['samtools'],
        'fastq',
        UNMAPPED_BAM,
        '-o', RESULT_UNMAPPED_FASTQ,
        '-0', RESULT_UNMAPPED_FASTQ
    ]

    cds_unmapped_samtools_view = [
        program_paths['samtools'],
        'view',
        '-h',
        '-b',
        '-f', '4',
        '-o', CDS_UNMAPPED_BAM,
        CDS_INTERMEDIATE_SAM
    ]

    cds_unmapped_bedtools_bamtofastq_command = [
        program_paths['samtools'],
        'fastq',
        CDS_UNMAPPED_BAM,
        '-o', CDS_UNMAPPED_FASTQ,
        '-0', CDS_UNMAPPED_FASTQ
    ]

    print('====> Building bowtie2 Index...')
    run(bowtie2_build_index, capture_output=not verbose)
    if cds is not None:
        run(cds_build_index, capture_output=not verbose)

    print('====> Aligning reads to the genome...')
    run(bowtie2_align_reads, capture_output=not verbose)

    run(samtools_view, capture_output=not verbose) 
    run(bedtools_bamtofastq_command, capture_output=not verbose)
    run(unmapped_samtools_view, capture_output=not verbose)
    run(unmapped_bedtools_bamtofastq_command, capture_output=not verbose)
    if cds is not None:
        run(cds_align_reads, capture_output=not verbose)

    print('====> Converting to FASTQ...')
    if cds is not None:
        run(cds_samtools_view, capture_output=not verbose)
        run(cds_bedtools_bamtofastq_command, capture_output=not verbose)
        run(cds_unmapped_samtools_view, capture_output=not verbose)
        run(cds_unmapped_bedtools_bamtofastq_command, capture_output=not verbose)

        copy2(CDS_UNMAPPED_FASTQ, FINAL_UNMAPPED_FASTQ)
    else:
        copy2(RESULT_UNMAPPED_FASTQ, FINAL_UNMAPPED_FASTQ)

    SeqIO.write(into_seqrecord(merge_genome_cds_hits(RESULT_FASTQ, CDS_RESULT_FASTQ)), FINAL_FASTQ, 'fastq')
    create_stats_table(small_rnas, output)

    if not keep_intermediate:
        shutil.rmtree(INDEX_DIRECTORY)

        os.remove(INTERMEDIATE_SAM)
        os.remove(INTERMEDIATE_BAM)
        os.remove(UNMAPPED_BAM)
        os.remove(RESULT_FASTQ)
        os.remove(RESULT_UNMAPPED_FASTQ)

        if cds is not None:
            os.remove(CDS_INTERMEDIATE_SAM)
            os.remove(CDS_INTERMEDIATE_BAM)
            os.remove(CDS_UNMAPPED_BAM)
            os.remove(CDS_UNMAPPED_FASTQ)
            os.remove(CDS_RESULT_FASTQ)

    return FINAL_FASTQ

def create_stats_table(smallrna, output_dir):
    '''
    Count sequences to create the statistics file
    '''

    input_reads = 0
    for read in SeqIO.parse(smallrna, 'fastq'):
        input_reads += 1

    overall_mapped = 0
    for read in SeqIO.parse(os.path.join(output_dir, 'mapped_sequences.fastq'), 'fastq'):
        overall_mapped += 1

    overall_unmapped = input_reads - overall_mapped

    genome_mapped = 0
    for read in SeqIO.parse(os.path.join(output_dir, 'genome_mapped_sequences.fastq'), 'fastq'):
        genome_mapped += 1

    cds_mapped = 0
    try:
        for read in SeqIO.parse(os.path.join(output_dir, 'cds_sequences.fastq'), 'fastq'):
            cds_mapped += 1
    except:
        pass

    cds_unmapped = input_reads - cds_mapped
    genome_unmapped = input_reads - genome_mapped

    with open(os.path.join(output_dir, 'alignment_report.tsv'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')

        writer.writerow(['Total Reads', input_reads])
        writer.writerow(['Total Mapped', overall_mapped])
        writer.writerow(['Percentage Mapped', overall_mapped / input_reads * 100])
        writer.writerow(['Total Unmapped', overall_unmapped])
        writer.writerow(['Percentage Unmapped', overall_unmapped / input_reads * 100])
        writer.writerow(['Mapped to Genome', genome_mapped])
        writer.writerow(['Percentage Mapped to Genome', genome_mapped / input_reads * 100])
        writer.writerow(['Unmapped to Genome', genome_unmapped])
        writer.writerow(['Percentage Unmapped to Genome', genome_unmapped / input_reads * 100])
        writer.writerow(['Mapped to CDS', cds_mapped])
        writer.writerow(['Percentage Mapped to CDS (but not to genome)', cds_mapped / input_reads * 100])
        writer.writerow(['Unmapped to CDS or genome', cds_unmapped])
        writer.writerow(['Percentage Unmapped to CDS or genome', cds_unmapped / input_reads * 100])

def bin_rna_size(rna_file, min_length, max_length, output, verbose=False):
    '''
    Bin the RNAs in the RNA file into new files by length and first base
    '''
    BIN_LENGTH_DIRECTORY = os.path.join(output, 'binned_length_rna')
    if not os.path.exists(BIN_LENGTH_DIRECTORY):
        os.makedirs(BIN_LENGTH_DIRECTORY)

    BIN_BASE_DIRECTORY = os.path.join(output, 'binned_group_rna')
    if not os.path.exists(BIN_BASE_DIRECTORY):
        os.makedirs(BIN_BASE_DIRECTORY)

    print('====> Sorting RNA into arrays by length and first base...')

    length_firstbase = defaultdict(lambda: defaultdict(lambda: []))
    sequence_counts = defaultdict(lambda: 0)

    if os.path.exists(os.path.join(output, 'unmapped_sequences.fastq')):
        unmapped = 0
        for seq in SeqIO.parse(os.path.join(output, 'unmapped_sequences.fastq'), 'fastq'):
            unmapped += 1

        sequence_counts['__not_aligned'] = unmapped
    else:
        sequence_counts['__not_aligned'] = 0

    for seq in SeqIO.parse(rna_file, 'fastq'):
        first_base = seq.seq[0].replace('T', 'U')
        if len(seq) >= min_length and len(seq) <= max_length:
            sequence_counts[str(seq.seq)] += 1
            length_firstbase[len(seq)][first_base].append(seq)

    for length in length_firstbase:
        length_rnas = [x for xs in length_firstbase[length].values() for x in xs]
        SeqIO.write(length_rnas, os.path.join(BIN_LENGTH_DIRECTORY, f'length{length}.fastq'), 'fastq')

        for first_base in length_firstbase[length]:
            SeqIO.write(length_firstbase[length][first_base], os.path.join(BIN_BASE_DIRECTORY, f'{length}{first_base}.fastq'), 'fastq')

    table_path = os.path.join(output, 'rna_length_report.csv')
    counts_path = os.path.join(output, 'counts.tab')
    with open(counts_path, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for sequence in sequence_counts:
            writer.writerow([sequence, sequence_counts[sequence]])

    with open(table_path, 'w') as table_file:
        writer = csv.writer(table_file)
        writer.writerow(['RNA Length', 'A', 'C', 'G', 'U', 'N', 'All Frequency'])

        for length in range(min_length, max_length + 1):
            length_counts = [
                len(length_firstbase[length]['A']), len(length_firstbase[length]['C']),
                len(length_firstbase[length]['G']), len(length_firstbase[length]['U']),
                len(length_firstbase[length]['N'])
            ]

            writer.writerow([length] + length_counts + [sum(length_counts)])

    return table_path

def graph_length(path_to_table, output):
    '''
    Create the graph showing the length and starting base of all the RNA
    '''
    with open(path_to_table) as f:
        reader = csv.reader(f)
        next(reader)

        labels = []
        bases = {'A': [], 'C': [], 'G': [], 'T': [], 'N': []}
        
        for line in reader:
            labels.append(str(line[0]))

            bases['A'].append(int(line[1]))
            bases['C'].append(int(line[2]))
            bases['G'].append(int(line[3]))
            bases['T'].append(int(line[4]))
            bases['N'].append(int(line[5]))

        total_rna = sum(bases['A']) + sum(bases['C']) + sum(bases['G']) + sum(bases['T']) + sum(bases['N'])
        for key in bases.keys():
            bases[key] = np.array(bases[key]) / total_rna * 100

        with open(os.path.join(output, 'baseplot_data.csv'), 'w') as f:
            writer = csv.writer(f)

            writer.writerow(['RNA Length'] + labels)
            writer.writerow(['A'] + list(bases['A']))
            writer.writerow(['T'] + list(bases['T']))
            writer.writerow(['C'] + list(bases['C']))
            writer.writerow(['G'] + list(bases['G']))
            writer.writerow(['N'] + list(bases['N']))

        fig, ax = plt.subplots(figsize=(11.69, 8.27))

        ax.bar(labels, bases['A'], label='A', align='edge', width=0.8)
        ax.bar(labels, bases['C'], bottom=bases['A'], label='C', align='edge', width=0.8)
        ax.bar(labels, bases['G'], bottom=bases['A'] + bases['C'], label='G', align='edge', width=0.8)
        ax.bar(labels, bases['T'], bottom=bases['A'] + bases['C'] + bases['G'], label='U', align='edge', width=0.8)
        ax.bar(labels, bases['N'], bottom=bases['A'] + bases['C'] + bases['G'] + bases['T'], label='N', align='edge', width=0.8)

        ax.set_xticks(list(range(0, len(labels))))
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_xlabel('Length of small RNA')
        ax.set_ylabel('Percentage of small RNA')

        plt.legend(title='First Nucliotide')

        plt.savefig(os.path.join(output, 'baseplot.png'))