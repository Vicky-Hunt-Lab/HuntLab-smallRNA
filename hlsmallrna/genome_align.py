import os
import os.path
import csv

from subprocess import run

import numpy as np

from Bio import SeqIO
from matplotlib import pyplot as plt

from .config import get_config_key, mkdir_if_not_exists, do_log

def align_to_genome(genome, small_rnas, quiet=0):
    '''
    Align the small RNAs to the genome and filter out any that are unsuccessful
    '''

    RESULT_FASTQ = os.path.join(get_config_key('general', 'output_directory'), 'mapped_sequences.fastq')
    INDEX_DIRECTORY = os.path.join(get_config_key('general', 'output_directory'), 'bbmap_index')

    mkdir_if_not_exists(INDEX_DIRECTORY)

    bbmap_build_index = [
        get_config_key('cli-tools', 'bbmap', 'path_to_bbmap'),
        'ref=' + genome,
        'path=' + INDEX_DIRECTORY
    ]

    bbmap_build_index = bbmap_build_index + get_config_key('cli-tools', 'bbmap', 'bbmap_index_params')

    bbmap_align_reads = [
        get_config_key('cli-tools', 'bbmap', 'path_to_bbmap'),
        'in=' + small_rnas,
        'out=' + RESULT_FASTQ,
        'path=' + INDEX_DIRECTORY
    ]

    bbmap_align_reads = bbmap_align_reads + get_config_key('cli-tools', 'bbmap', 'bbmap_align_params')

    if get_config_key('cli-tools', 'bbmap', 'bbmap_pass_threads'):
        threads = get_config_key('general', 'threads')

        bbmap_build_index.append('threads=' + str(threads))
        bbmap_align_reads.append('threads=' + str(threads))

    do_log(quiet, '====> Building BBMap Index')
    run(bbmap_build_index, capture_output=(quiet != 0))
    do_log(quiet, '====> Aligning reads to the genome')
    run(bbmap_align_reads, capture_output=(quiet != 0))

    return RESULT_FASTQ

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