import os
import os.path
import csv

from subprocess import run

import pysam

from Bio import SeqIO

from config import get_config_key, mkdir_if_not_exists

def revcomp_input_file(smallRNA):
    '''
    Create a file containing the reverse complement of a file of small RNA
    '''

    print('====> Reverse complimenting RNA...')
    PATH_TO_REVCOMP = os.path.join(get_config_key('general', 'output_directory'), 'revcomp_rna.fastq')

    def do_revcomp(seq):
        revcomp = seq.reverse_complement()
        revcomp.id = seq.id

        return revcomp

    seqs = SeqIO.parse(smallRNA, 'fastq')
    revcomps = map(do_revcomp, seqs)
    SeqIO.write(revcomps, PATH_TO_REVCOMP, 'fastq')

    return PATH_TO_REVCOMP


def find_targets(smallRNA, possible_target_list, min_seq_length=2, quiet=False):
    '''
    Align the small RNA against the lists of possible targets with bowtie2 and
    analyse the output
    '''

    CWD = os.getcwd()
    INDEX_DIR = os.path.join(get_config_key('general', 'output_directory'), 'bowtie_indexes')
    SAM_DIR = os.path.join(CWD, get_config_key('general', 'output_directory'), 'target_alignments')

    sam_files = []

    mkdir_if_not_exists(SAM_DIR)
    mkdir_if_not_exists(INDEX_DIR)
    os.chdir(INDEX_DIR)

    for target in possible_target_list:
        print(f'====> Builing index for {target}')
        INDEX_NAME = os.path.basename(target) + '_index'

        bowtie2_build_command = [
            get_config_key('cli-tools', 'bowtie2', 'path_to_bowtie2_build'),
            os.path.join(CWD, target),
            INDEX_NAME
        ]

        bowtie2_build_command = bowtie2_build_command + get_config_key('cli-tools', 'bowtie2', 'bowtie2_build_params')

        run(bowtie2_build_command, capture_output=quiet)

        print(f'====> Aligning small RNA against {target}')

        sam_files.append(os.path.join(SAM_DIR, INDEX_NAME + '.sam'))
        bowtie2_align_command = [
            get_config_key('cli-tools', 'bowtie2', 'path_to_bowtie2'),
            '-L', str(min_seq_length),
            '--no-1mm-upfront',
            '--score-min', 'L,0,0',
            '--end-to-end',
            '--norc',
            '-M', '0',
            '-x', INDEX_NAME,
            '-U', os.path.join(CWD, smallRNA),
            '-S', sam_files[-1]
        ]

        bowtie2_align_command = bowtie2_align_command + get_config_key('cli-tools', 'bowtie2', 'bowtie2_params')

        run(bowtie2_align_command, capture_output=quiet)

    os.chdir(CWD)

    return sam_files

def build_summery_files(sam_files):
    '''
    Takes the alignment SAM and extracts the sequences into a summary file and
    a set of FASTA files
    '''
    
    with open(os.path.join(get_config_key('general', 'output_directory'), 'rna_target_list.csv'), 'w') as tablefile:
        writer = csv.writer(tablefile)
        writer.writerow(['Query Sequence', 'Target File', 'Target Sequence'])

        for filename in sam_files:
            samfile = pysam.AlignmentFile(filename, 'r')
            fileheader = samfile.header

            remove_duplicate_command = [
                get_config_key('cli-tools', 'samtools', 'path_to_samtools'),
                'view',
                '-o', filename + '.nodup.sam',
                filename
            ]

            remove_duplicate_command = remove_duplicate_command + get_config_key('cli-tools', 'samtools', 'samtools_view_params')

            run(remove_duplicate_command)
            filename = filename + '.nodup.sam'

            for read in samfile:
                if not read.is_unmapped:
                    query_name = read.query_name
                    r_name = fileheader.get_reference_name(read.reference_id)

                    writer.writerow([query_name, os.path.abspath(filename), r_name])

            bam_to_fastq_command = [
                get_config_key('cli-tools', 'samtools', 'path_to_samtools'),
                'fastq',
                filename
            ]

            bam_to_fastq_command = bam_to_fastq_command + get_config_key('cli-tools', 'samtools', 'samtools_fastq_params')

            result = run(bam_to_fastq_command, capture_output=True)
            with open(filename + '.fastq', 'wb') as f:
                f.write(result.stdout)