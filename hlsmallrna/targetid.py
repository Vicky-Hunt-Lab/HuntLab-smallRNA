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
import os
import os.path
import csv

from subprocess import run
from collections import defaultdict

import pysam

from Bio import SeqIO
from Bio.Seq import Seq

def check_ids_unique(path_to_file, format):
    '''
    Make sure that the sequences in a FASTA or FASTQ file are unique
    '''
    seen_ids = set()
    for seq in SeqIO.parse(path_to_file, format):
        if seq.id in seen_ids:
            return False
        else:
            seen_ids.add(seq.id)

    return True

def seqio_revcomp_multiple(fastq_files):
    '''
    Merge SeqIO.parse over multiple files and reverse complement sequences
    '''
    seen_ids = set()

    for filepath in fastq_files:
        for seq in SeqIO.parse(filepath, 'fastq'):
            if seq.id in seen_ids:
                print(f'Error: duplicate ID {seq.id} found in your small RNA files, make IDs unique and try again')
                exit(1)
            else:
                seen_ids.add(seq.id)

            revcomp = seq.reverse_complement()
            revcomp.id = seq.id

            yield revcomp

def merge_and_revcomp_input_files(smallRNA_files, output):
    '''
    Create a file containing the reverse complement of a file of small RNA
    '''

    print('====> Reverse complimenting RNA...')
    PATH_TO_REVCOMP = os.path.join(output, 'revcomp_rna.fastq')

    SeqIO.write(seqio_revcomp_multiple(smallRNA_files), PATH_TO_REVCOMP, 'fastq')

    return PATH_TO_REVCOMP


def find_targets(smallRNA, possible_target_list, output, program_paths, threads=4, min_seq_length=2, mismatches_allowed=0, verbose=False):
    '''
    Align the small RNA against the lists of possible targets with bowtie2 and
    analyse the output
    '''

    CWD = os.getcwd()
    INDEX_DIR = os.path.join(output, 'bowtie_indexes')
    SAM_DIR = os.path.join(CWD, output, 'target_alignments')

    sam_files = []

    if not os.path.exists(SAM_DIR):
        os.makedirs(SAM_DIR)

    if not os.path.exists(INDEX_DIR):
        os.makedirs(INDEX_DIR)

    os.chdir(INDEX_DIR)

    for target in possible_target_list:
        if not check_ids_unique(os.path.join(CWD, target), 'fasta'):
            print(f'Error: duplicate ID found in your target file: {target}. Make IDs unique and try again')
            exit(1)

        print(f'====> Builing index for {target}...')
        INDEX_NAME = os.path.basename(target) + '_index'

        bowtie2_build_command = [
            program_paths['bowtie2-build'],
            '--threads', str(threads),
            os.path.join(CWD, target),
            INDEX_NAME
        ]

        run(bowtie2_build_command, capture_output=not verbose)

        print(f'====> Aligning small RNA against {target}...')

        sam_files.append(os.path.join(SAM_DIR, INDEX_NAME + '.sam'))

        if mismatches_allowed > 0:
            bowtie2_align_command = [
                program_paths['bowtie2'],
                '--threads', str(threads),
                '-L', str(min_seq_length),
                '--score-min', 'L,-' + str(mismatches_allowed) + ',0',
                '--end-to-end',
                '--norc',
                '--mp', '1,1',
                '--ignore-quals',
                '--rdg', '100,100',
                '--rfg', '100,100',
                '-x', INDEX_NAME,
                '-U', os.path.join(CWD, smallRNA),
                '-S', sam_files[-1],
                '-a'
            ]
        else:
            bowtie2_align_command = [
                program_paths['bowtie2'],
                '--threads', str(threads),
                '-L', str(min_seq_length),
                '--no-1mm-upfront',
                '--score-min', 'L,0,0',
                '--end-to-end',
                '--norc',
                '-M', '0',
                '-x', INDEX_NAME,
                '-U', os.path.join(CWD, smallRNA),
                '-S', sam_files[-1],
                '-a'
            ]

        run(bowtie2_align_command, capture_output=not verbose)

    os.chdir(CWD)

    return sam_files

def build_summary_files(sam_files, output, program_paths, verbose=False):
    '''
    Takes the alignment SAM and extracts the sequences into a summary file and
    a set of FASTA files
    '''
    print('====> Removing unaligend sequences and building summary files...')

    target_rows = defaultdict(lambda: [])
    for filename in sam_files:
            samfile = pysam.AlignmentFile(filename, 'r')
            fileheader = samfile.header

            for read in samfile:
                if not read.is_unmapped:
                    query_name = read.query_name
                    query_sequence = Seq(read.query_sequence).reverse_complement()
                    r_name = fileheader.get_reference_name(read.reference_id)
                    start = read.reference_start
                    end = read.reference_end

                    if read.is_reverse:
                        strand = 'antisense'
                    else:
                        strand = 'sense'

                    length = len(query_sequence)
                    first_base = query_sequence[0].replace('T', 'U')
                    target_rows[f'{length}{first_base}'].append([query_name, query_sequence, os.path.abspath(filename), r_name, start, end, strand])


    TRAGETS_DIR = os.path.join(output, 'targets_by_group')

    if not os.path.exists(TRAGETS_DIR):
        os.makedirs(TRAGETS_DIR)

    for target in target_rows:
        with open(os.path.join(TRAGETS_DIR, f'{target}_targets.tsv'), 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['Query Name', 'Query Sequence', 'Target File', 'Target Sequence', 'Start Coordinate', 'End Coordinate', 'Strand'])
            writer.writerows(target_rows[target])
    
    with open(os.path.join(output, 'rna_target_list.tsv'), 'w') as tablefile:
        writer = csv.writer(tablefile, delimiter='\t')
        writer.writerow(['Query Name', 'Query Sequence', 'Target File', 'Target Sequence', 'Start Coordinate', 'End Coordinate', 'Strand'])
        for target in target_rows:
            writer.writerows(target_rows[target])

    try:
        remove_duplicate_command = [
            program_paths['samtools'],
            'view',
            '-h', 
            '-F', '256',
            '-F', '4',
            '-o', filename + '.nodup.sam',
            filename
        ]
    
        result = run(remove_duplicate_command, capture_output=not verbose)
        result.check_returncode()
        filename = filename + '.nodup.sam'

        bam_to_fastq_command = [
            program_paths['samtools'],
            'fastq',
            filename
        ]

        result = run(bam_to_fastq_command, capture_output=True)
        result.check_returncode()
        with open(filename + '.fastq', 'wb') as f:
            f.write(result.stdout)

    except:
        print('====> SAMTOOLS convertion failed, skipping creating filtered and FASTQ files...')