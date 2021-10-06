import os
import os.path

from subprocess import run
from collections import defaultdict

from Bio import SeqIO

from config import get_config_key

def align_to_genome(genome, small_rnas, quiet=False):
    '''
    Align the small RNAs to the genome and filter out any that are unsuccessful
    '''

    CWD = os.getcwd()
    INDEX_DIR = os.path.join(CWD, get_config_key('general', 'output_directory'), 'bowtie-indexs')
    BOWTIE2_DIR = os.path.join(CWD, get_config_key('general', 'output_directory'), 'align_and_filter')
    ALIGNMENT_SAM = os.path.join(BOWTIE2_DIR, 'aligned_rnas.sam')
    VIEW_BAM = os.path.join(BOWTIE2_DIR, 'view_rnas.bam')
    SORTED_BAM = os.path.join(BOWTIE2_DIR, 'sorted_rnas.bam')
    RESULT_FASTQ = os.path.join(BOWTIE2_DIR, 'filtered_result.fastq')

    genome_path = os.path.join(CWD, genome)
    small_rna_path = os.path.join(CWD, small_rnas)
 
    os.mkdir(INDEX_DIR)
    os.mkdir(BOWTIE2_DIR)
    os.chdir(INDEX_DIR)

    bowtie_build_command = [
        get_config_key('cli-tools', 'bowtie2', 'path_to_bowtie2_build'),
        genome_path,
        'genome_index'
    ]

    bowtie_build_command = bowtie_build_command + get_config_key('cli-tools', 'bowtie2', 'bowtie2_build_params')

    # TODO: make bowtie only do one alignment
    
    bowtie_align_command = [
        get_config_key('cli-tools', 'bowtie2', 'path_to_bowtie2'),
        '-x', 'genome_index',
        '-k', '1',
        '-U', small_rna_path,
        '-S', ALIGNMENT_SAM
    ]

    bowtie_align_command = bowtie_align_command + get_config_key('cli-tools', 'bowtie2', 'bowtie2_params')

    samtools_view_command = [
        get_config_key('cli-tools', 'samtools', 'path_to_samtools'),
        'view',
        ALIGNMENT_SAM,
        '-F', '2',
        '-o', VIEW_BAM
    ]

    samtools_view_command = samtools_view_command + get_config_key('cli-tools', 'samtools', 'samtools_view_params')

    samtools_sort_command = [
        get_config_key('cli-tools', 'samtools', 'path_to_samtools'),
        'sort',
        '-n', VIEW_BAM,
        '-o', SORTED_BAM
    ]

    samtools_sort_command = samtools_sort_command + get_config_key('cli-tools', 'samtools', 'samtools_sort_params')

    bamToFastq_command = [
        get_config_key('cli-tools', 'bedtools', 'path_to_bedtools'),
        'bamtofastq',
        '-i', SORTED_BAM,
        '-fq', RESULT_FASTQ
    ]

    bamToFastq_command = bamToFastq_command + get_config_key('cli-tools', 'bedtools', 'bedtools_bamToFastq_params')

    print('====> Running step 1/4, Building bowtie index')
    result = run(bowtie_build_command, capture_output=quiet)
    result.check_returncode()

    print('====> Running step 2/4, Aligning RNA to genome')
    result = run(bowtie_align_command, capture_output=quiet)
    result.check_returncode()

    print('====> Running step 3/4, Sorting alignments')
    result = run(samtools_view_command, capture_output=quiet)
    result.check_returncode()

    result = run(samtools_sort_command, capture_output=quiet)
    result.check_returncode()

    print('====> Running step 4/4, Converting aligned RNA to FASTQ')
    result = run(bamToFastq_command, capture_output=quiet)
    result.check_returncode()

    os.chdir(CWD)

    return RESULT_FASTQ

def bin_rna_size(rna_file):
    '''
    Bin the RNAs in the RNA file into new files by length
    '''

    BINS_DIRECTORY = os.path.join(get_config_key('general', 'output_directory'), 'binned_rna')
    os.mkdir(BINS_DIRECTORY)

    print('====> Sorting RNA into arrays by length')
    files = defaultdict(lambda: [])
    rnas = SeqIO.parse(rna_file, 'fastq')
    for rna_seq in rnas:
        files[len(rna_seq)].append(rna_seq)

    print('====> Writing RNA to new files')
    for key in files.keys():
        filename = os.path.join(BINS_DIRECTORY, 'length' + str(key) + '.fastq')
        SeqIO.write(files[key], filename, 'fastq')
