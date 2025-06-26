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
from shutil import copy2

import os
import os.path
import csv

from subprocess import run

import numpy as np

from Bio import SeqIO, SeqRecord, Seq
from matplotlib import pyplot as plt

from .config import get_config_key, mkdir_if_not_exists, do_log

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

def remove_symbols_from_header(fasta):
    '''
    Remove @s from the FASTA header before doing the alignment
    '''
    for seq in SeqIO.parse(fasta, 'fasta'):
        seq.id = seq.id.replace('@', '_')

        yield seq

def align_to_genome(genome, small_rnas, cds, quiet=0, threads=4, small_rna_filetype='fastq', mismatches=0):
    '''
    Align the small RNAs to the genome and filter out any that are unsuccessful
    '''

    INTERMEDIATE_SAM = os.path.join(get_config_key('general', 'output_directory'), 'mapped_sequences.sam')
    INTERMEDIATE_BAM = os.path.join(get_config_key('general', 'output_directory'), 'mapped_sequences.bam')
    UNMAPPED_BAM = os.path.join(get_config_key('general', 'output_directory'), 'unmapped_sequences.bam')
    RESULT_FASTQ = os.path.join(get_config_key('general', 'output_directory'), 'genome_mapped_sequences.fastq')
    RESULT_UNMAPPED_FASTQ = os.path.join(get_config_key('general', 'output_directory'), 'genome_unmapped_sequences.fastq')
    INDEX_DIRECTORY = os.path.join(get_config_key('general', 'output_directory'), 'bbmap_index')

    CDS_INTERMEDIATE_SAM = os.path.join(get_config_key('general', 'output_directory'), 'cds_sequences.sam')
    CDS_INTERMEDIATE_BAM = os.path.join(get_config_key('general', 'output_directory'), 'cds_sequences.bam')
    CDS_UNMAPPED_BAM = os.path.join(get_config_key('general', 'output_directory'), 'cds_unmapped_sequences.bam')
    CDS_UNMAPPED_FASTQ = os.path.join(get_config_key('general', 'output_directory'), 'cds_unmapped_sequences.fastq')
    CDS_RESULT_FASTQ = os.path.join(get_config_key('general', 'output_directory'), 'cds_sequences.fastq')

    FINAL_FASTQ = os.path.join(get_config_key('general', 'output_directory'), 'mapped_sequences.fastq')
    FINAL_UNMAPPED_FASTQ = os.path.join(get_config_key('general', 'output_directory'), 'unmapped_sequences.fastq')

    mkdir_if_not_exists(INDEX_DIRECTORY)

    if small_rna_filetype == 'fasta':
        NEW_SMALL_RNAS = os.path.join(get_config_key('general', 'output_directory'), 'corrected_headers.fasta')
        SeqIO.write(remove_symbols_from_header(small_rnas), NEW_SMALL_RNAS, 'fasta')

        small_rnas = NEW_SMALL_RNAS

    bbmap_build_index = [
        get_config_key('cli-tools', 'bowtie2', 'path_to_bowtie2_build'),
        '--threads', str(threads),
        genome,
        os.path.join(INDEX_DIRECTORY, 'genome_index')
    ]

    bbmap_build_index = bbmap_build_index + get_config_key('cli-tools', 'bowtie2', 'bowtie2_build_params')

    cds_build_index = [
        get_config_key('cli-tools', 'bowtie2', 'path_to_bowtie2_build'),
        '--threads', str(threads),
        cds,
        os.path.join(INDEX_DIRECTORY, 'cds_index')
    ]

    cds_build_index = cds_build_index + get_config_key('cli-tools', 'bowtie2', 'bowtie2_build_params')

    if mismatches == 0:
        bbmap_align_reads = [
                get_config_key('cli-tools', 'bowtie2', 'path_to_bowtie2'),
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
        bbmap_align_reads = [
            get_config_key('cli-tools', 'bowtie2', 'path_to_bowtie2'),
            '--threads', str(threads),
            '-L', '18',
            '--score-min', 'L,' + str(-(mismatches - 1)) + ',0',
            '--end-to-end',
            '--mp', '1,1',
            '--ignore-quals',
            '--rdg', '9,1',
            '--rfg', '9,1',
            '-x', os.path.join(INDEX_DIRECTORY, 'genome_index'),
            '-U', small_rnas,
            '-S', INTERMEDIATE_SAM
        ]

    if small_rna_filetype == 'fasta':
        bbmap_align_reads.append('-f')

    bbmap_align_reads = bbmap_align_reads + get_config_key('cli-tools', 'bowtie2', 'bowtie2_params')


    if mismatches == 0:
        cds_align_reads = [
                get_config_key('cli-tools', 'bowtie2', 'path_to_bowtie2'),
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
            get_config_key('cli-tools', 'bowtie2', 'path_to_bowtie2'),
            '--threads', str(threads),
            '-L', '18',
            '--score-min', 'L,' + str(-(mismatches - 1)) + ',0',
            '--end-to-end',
            '--mp', '1,1',
            '--ignore-quals',
            '--rdg', '9,1',
            '--rfg', '9,1',
            '-x', os.path.join(INDEX_DIRECTORY, 'cds_index'),
            '-U', RESULT_UNMAPPED_FASTQ,
            '-S', CDS_INTERMEDIATE_SAM
        ]

    cds_align_reads = cds_align_reads + get_config_key('cli-tools', 'bowtie2', 'bowtie2_params')

    samtools_view = [
        get_config_key('cli-tools', 'samtools', 'path_to_samtools'),
        'view',
        '-h',
        '-b',
        '-F', '4',
        '-o', INTERMEDIATE_BAM,
        INTERMEDIATE_SAM
    ]

    samtools_view = samtools_view + get_config_key('cli-tools', 'samtools', 'samtools_view_params')

    bedtools_bamtofastq_command = [
        get_config_key('cli-tools', 'samtools', 'path_to_samtools'),
        'fastq',
        INTERMEDIATE_BAM,
        '-o', RESULT_FASTQ,
        '-0', RESULT_FASTQ
    ]

    cds_samtools_view = [
        get_config_key('cli-tools', 'samtools', 'path_to_samtools'),
        'view',
        '-h',
        '-b',
        '-F', '4',
        '-o', CDS_INTERMEDIATE_BAM,
        CDS_INTERMEDIATE_SAM
    ]

    cds_samtools_view = cds_samtools_view + get_config_key('cli-tools', 'samtools', 'samtools_view_params')

    cds_bedtools_bamtofastq_command = [
        get_config_key('cli-tools', 'samtools', 'path_to_samtools'),
        'fastq',
        CDS_INTERMEDIATE_BAM,
        '-o', CDS_RESULT_FASTQ,
        '-0', CDS_RESULT_FASTQ
    ]
    
    unmapped_samtools_view = [
        get_config_key('cli-tools', 'samtools', 'path_to_samtools'),
        'view',
        '-h',
        '-b',
        '-f', '4',
        '-o', UNMAPPED_BAM,
        INTERMEDIATE_SAM
    ]

    unmapped_bedtools_bamtofastq_command = [
        get_config_key('cli-tools', 'samtools', 'path_to_samtools'),
        'fastq',
        UNMAPPED_BAM,
        '-o', RESULT_UNMAPPED_FASTQ,
        '-0', RESULT_UNMAPPED_FASTQ
    ]

    cds_unmapped_samtools_view = [
        get_config_key('cli-tools', 'samtools', 'path_to_samtools'),
        'view',
        '-h',
        '-b',
        '-f', '4',
        '-o', CDS_UNMAPPED_BAM,
        CDS_INTERMEDIATE_SAM
    ]

    cds_unmapped_bedtools_bamtofastq_command = [
        get_config_key('cli-tools', 'samtools', 'path_to_samtools'),
        'fastq',
        CDS_UNMAPPED_BAM,
        '-o', CDS_UNMAPPED_FASTQ,
        '-0', CDS_UNMAPPED_FASTQ
    ]

    do_log(quiet, '====> Building BBMap Index')
    run(bbmap_build_index, capture_output=(quiet != 0))
    if cds is not None:
        run(cds_build_index, capture_output=(quiet != 0))

    do_log(quiet, '====> Aligning reads to the genome')
    run(bbmap_align_reads, capture_output=(quiet != 0))

    run(samtools_view, capture_output=(quiet != 0)) 
    run(bedtools_bamtofastq_command, capture_output=(quiet != 0))
    run(unmapped_samtools_view, capture_output=(quiet != 0))
    run(unmapped_bedtools_bamtofastq_command, capture_output=(quiet != 0))
    if cds is not None:
        run(cds_align_reads, capture_output=(quiet != 0))

    do_log(quiet, '====> Converting to FASTQ')
    if cds is not None:
        run(cds_samtools_view, capture_output=(quiet != 0))
        run(cds_bedtools_bamtofastq_command, capture_output=(quiet != 0))
        run(cds_unmapped_samtools_view, capture_output=(quiet != 0))
        run(cds_unmapped_bedtools_bamtofastq_command, capture_output=(quiet != 0))

        copy2(CDS_UNMAPPED_FASTQ, FINAL_UNMAPPED_FASTQ)
    else:
        copy2(RESULT_UNMAPPED_FASTQ, FINAL_UNMAPPED_FASTQ)

    SeqIO.write(into_seqrecord(merge_genome_cds_hits(RESULT_FASTQ, CDS_RESULT_FASTQ)), FINAL_FASTQ, 'fastq')
    create_stats_table(small_rnas, get_config_key('general', 'output_directory'), small_rna_filetype=small_rna_filetype)

    return FINAL_FASTQ

def create_stats_table(smallrna, output_dir, small_rna_filetype='fastq'):
    '''
    Count sequences to create the statistics file
    '''

    input_reads = 0
    for read in SeqIO.parse(smallrna, small_rna_filetype):
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

    with open(os.path.join(output_dir, 'report.tsv'), 'w') as f:
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
        writer.writerow(['Percentage Mapped to CDS', cds_mapped / input_reads * 100])
        writer.writerow(['Unmapped to CDS', cds_unmapped])
        writer.writerow(['Percentage Unmapped to CDS', cds_unmapped / input_reads * 100])

def bin_rna_size(rna_file, min_length, max_length, quiet=0):
    '''
    Bin the RNAs in the RNA file into new files by length
    '''

    BINS_DIRECTORY = os.path.join(get_config_key('general', 'output_directory'), 'binned_rna')
    mkdir_if_not_exists(BINS_DIRECTORY)

    do_log(quiet, '====> Sorting RNA into arrays by length')
    rnas = sorted(SeqIO.parse(rna_file, 'fastq'), key=lambda x: len(x))

    table_path = os.path.join(get_config_key('general', 'output_directory'), 'rna_length_report.csv')

    last_length = 0
    current_rnas = []
    start_base_count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    with open(table_path, 'w') as table_file:
        writer = csv.writer(table_file)
        writer.writerow(['RNA Length', 'A', 'C', 'G', 'U', 'All Frequency'])

        for rna_seq in rnas:
            if len(rna_seq) != last_length:
                if len(current_rnas) > 0 and last_length >= min_length and last_length <= max_length:
                    filename = os.path.join(BINS_DIRECTORY, 'length' + str(last_length) + '.fastq')
                    with open(filename, 'a') as f:
                        SeqIO.write(current_rnas, f, 'fastq')

                    writer.writerow([
                        last_length, 
                        start_base_count['A'], start_base_count['C'],
                        start_base_count['G'], start_base_count['T'],
                        len(current_rnas)
                    ])

                current_rnas = []
                start_base_count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

                last_length = len(rna_seq)

            current_rnas.append(rna_seq)
            try:
                start_base_count[rna_seq[0]] += 1
            except KeyError:
                do_log(quiet, f'Warning : An RNA started with ambiguous nt {rna_seq[0]}')

    return table_path

def graph_length(path_to_table):
    '''
    Create the graph showing the length and starting base of all the RNA
    '''
    
    with open(path_to_table) as f:
        reader = csv.reader(f)
        next(reader)

        labels = []
        bases = {'A': [], 'C': [], 'G': [], 'T': []}
        
        for line in reader:
            labels.append(str(line[0]))

            bases['A'].append(int(line[1]))
            bases['C'].append(int(line[2]))
            bases['G'].append(int(line[3]))
            bases['T'].append(int(line[4]))

        total_rna = sum(bases['A']) + sum(bases['C']) + sum(bases['G']) + sum(bases['T'])
        for key in bases.keys():
            bases[key] = np.array(bases[key]) / total_rna * 100

        with open(os.path.join(get_config_key('general', 'output_directory'), 'baseplot_data.csv'), 'w') as f:
            writer = csv.writer(f)

            writer.writerow(['RNA Length'] + labels)
            writer.writerow(['A'] + list(bases['A']))
            writer.writerow(['T'] + list(bases['T']))
            writer.writerow(['C'] + list(bases['C']))
            writer.writerow(['G'] + list(bases['G']))

        fig, ax = plt.subplots(figsize=(11.69, 8.27))

        ax.bar(labels, bases['A'], label='A', align='edge', width=0.8)
        ax.bar(labels, bases['C'], bottom=bases['A'], label='C', align='edge', width=0.8)
        ax.bar(labels, bases['G'], bottom=bases['A'] + bases['C'], label='G', align='edge', width=0.8)
        ax.bar(labels, bases['T'], bottom=bases['A'] + bases['C'] + bases['G'], label='U', align='edge', width=0.8)

        ax.set_xticklabels(labels, fontsize=7)
        ax.set_xlabel('Length of small RNA')
        ax.set_ylabel('Percentage of small RNA')

        plt.legend(title='First Nucliotide')

        plt.savefig(os.path.join(get_config_key('general', 'output_directory'), 'baseplot.png'))