import os
import os.path

from subprocess import run
from collections import defaultdict

from Bio import SeqIO

from config import get_config_key, mkdir_if_not_exists

def align_to_genome(genome, small_rnas, quiet=False):
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

    print('====> Building BBMap Index')
    run(bbmap_build_index, capture_output=quiet)
    print('====> Aligning reads to the genome')
    run(bbmap_align_reads, capture_output=quiet)

    return RESULT_FASTQ

def bin_rna_size(rna_file):
    '''
    Bin the RNAs in the RNA file into new files by length
    '''

    BINS_DIRECTORY = os.path.join(get_config_key('general', 'output_directory'), 'binned_rna')
    mkdir_if_not_exists(BINS_DIRECTORY)

    print('====> Sorting RNA into arrays by length')
    files = defaultdict(lambda: [])
    rnas = SeqIO.parse(rna_file, 'fastq')
    for rna_seq in rnas:
        files[len(rna_seq)].append(rna_seq)

    print('====> Writing RNA to new files')
    for key in files.keys():
        filename = os.path.join(BINS_DIRECTORY, 'length' + str(key) + '.fastq')
        SeqIO.write(files[key], filename, 'fastq')
