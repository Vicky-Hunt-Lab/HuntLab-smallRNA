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
import os
import os.path
import csv
import shutil

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

def make_fastqs_unique(fastq1, fastq2, output):
    '''
    Take 2 FastQ files and make the sequence in them unique
    '''
    unique_sequences = set()

    for sequence in SeqIO.parse(fastq1, 'fastq'):
        unique_sequences.add(str(sequence.seq))

    try:
        for sequence in SeqIO.parse(fastq2, 'fastq'):
            unique_sequences.add(str(sequence.seq))
    except FileNotFoundError:
        pass

    SeqIO.write(into_seqrecord(unique_sequences), output, 'fastq')

def align_to_genome(genome, small_rnas, cds, quiet=0):
    '''
    Align the small RNAs to the genome and filter out any that are unsuccessful
    '''

    INTERMEDIATE_SAM = os.path.join(get_config_key('general', 'output_directory'), 'mapped_sequences.sam')
    INTERMEDIATE_BAM = os.path.join(get_config_key('general', 'output_directory'), 'mapped_sequences.bam')
    RESULT_FASTQ = os.path.join(get_config_key('general', 'output_directory'), 'genome_mapped_sequences.fastq')
    INDEX_DIRECTORY = os.path.join(get_config_key('general', 'output_directory'), 'bbmap_index')

    CDS_INTERMEDIATE_SAM = os.path.join(get_config_key('general', 'output_directory'), 'cds_sequences.sam')
    CDS_INTERMEDIATE_BAM = os.path.join(get_config_key('general', 'output_directory'), 'cds_sequences.bam')
    CDS_RESULT_FASTQ = os.path.join(get_config_key('general', 'output_directory'), 'cds_sequences.fastq')

    FINAL_FASTQ = os.path.join(get_config_key('general', 'output_directory'), 'mapped_sequences.fastq')

    mkdir_if_not_exists(INDEX_DIRECTORY)

    bbmap_build_index = [
        get_config_key('cli-tools', 'bowtie2', 'path_to_bowtie2_build'),
        genome,
        os.path.join(INDEX_DIRECTORY, 'genome_index')
    ]

    bbmap_build_index = bbmap_build_index + get_config_key('cli-tools', 'bowtie2', 'bowtie2_build_params')

    cds_build_index = [
        get_config_key('cli-tools', 'bowtie2', 'path_to_bowtie2_build'),
        cds,
        os.path.join(INDEX_DIRECTORY, 'cds_index')
    ]

    cds_build_index = cds_build_index + get_config_key('cli-tools', 'bowtie2', 'bowtie2_build_params')

    bbmap_align_reads = [
                get_config_key('cli-tools', 'bowtie2', 'path_to_bowtie2'),
                '-L', '18',
                '-x', os.path.join(INDEX_DIRECTORY, 'genome_index'),
                '-U', small_rnas,
                '-S', INTERMEDIATE_SAM
    ]

    bbmap_align_reads = bbmap_align_reads + get_config_key('cli-tools', 'bowtie2', 'bowtie2_params')

    cds_align_reads = [
                get_config_key('cli-tools', 'bowtie2', 'path_to_bowtie2'),
                '-L', '18',
                '-x', os.path.join(INDEX_DIRECTORY, 'cds_index'),
                '-U', small_rnas,
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
        get_config_key('cli-tools', 'bedtools', 'path_to_bedtools'),
        'bamtofastq',
        '-i', INTERMEDIATE_BAM,
        '-fq', RESULT_FASTQ
    ]

    bedtools_bamtofastq_command = bedtools_bamtofastq_command + get_config_key('cli-tools', 'bedtools', 'bedtools_bamtofastq_params')

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
        get_config_key('cli-tools', 'bedtools', 'path_to_bedtools'),
        'bamtofastq',
        '-i', CDS_INTERMEDIATE_BAM,
        '-fq', CDS_RESULT_FASTQ
    ]

    cds_bedtools_bamtofastq_command = cds_bedtools_bamtofastq_command + get_config_key('cli-tools', 'bedtools', 'bedtools_bamtofastq_params')

    if get_config_key('cli-tools', 'bowtie2', 'bowtie2_pass_threads'):
        threads = get_config_key('general', 'threads')

        bbmap_build_index.append('--threads')
        bbmap_build_index.append(str(threads))
        bbmap_align_reads.append('--threads')
        bbmap_align_reads.append(str(threads))

        cds_build_index.append('--threads')
        cds_build_index.append(str(threads))
        cds_align_reads.append('--threads')
        cds_align_reads.append(str(threads))

    do_log(quiet, '====> Building BBMap Index')
    run(bbmap_build_index, capture_output=(quiet != 0))
    if cds is not None:
        run(cds_build_index, capture_output=(quiet != 0))

    do_log(quiet, '====> Aligning reads to the genome')
    run(bbmap_align_reads, capture_output=(quiet != 0))
    if cds is not None:
        run(cds_align_reads, capture_output=(quiet != 0))

    do_log(quiet, '====> Converting to FASTQ')
    run(samtools_view, capture_output=(quiet != 0)) 
    run(bedtools_bamtofastq_command, capture_output=(quiet != 0))
    if cds is not None:
        run(cds_samtools_view, capture_output=(quiet != 0))
        run(cds_bedtools_bamtofastq_command, capture_output=(quiet != 0))

    make_fastqs_unique(RESULT_FASTQ, CDS_RESULT_FASTQ, FINAL_FASTQ)

    return FINAL_FASTQ

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