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
import gzip
import shutil

import yaml

from Bio import SeqIO

KNOWN_ADAPTERS = {
    'qiagen': {
        '5': 'GTTCAGAGTTCTACAGTCCGACGATC',
        '3': 'AACTGTAGGCACCATCAAT'
    },
    'qiagen-rev': {
        '5': 'AACTGTAGGCACCATCAAT',
        '3': 'GTTCAGAGTTCTACAGTCCGACGATC'
    },
    'nebnext': {
        '5': None,
        '3': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    },
    'nebnext-rev': {
        '5': None,
        '3': 'GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT'
    }
}

def unzip_if_needed(filepath, output_dir, tmpfilename='temp.fq'):
    '''
    Check if filepath is gzipped
    '''
    with open(filepath, 'rb') as test_f:
        is_gzipped = test_f.read(2) == b'\x1f\x8b'

    if not is_gzipped:
        return filepath
    else:
        with gzip.open(filepath, 'rb') as file_in:
            with open(os.path.join(output_dir, tmpfilename), 'wb') as file_out:
                shutil.copyfileobj(file_in, file_out)

        return os.path.join(output_dir, tmpfilename)


def fasta_to_fastq_record(seqs):
    '''
    remove @s from the sequence IDs of a biopython object, replace with _ and
    add phread scores of 40
    '''
    for seq in seqs:
        if '@' in seq.id:
            seq.id = seq.id.replace('@', '_')

        seq.letter_annotations["phred_quality"] = [40] * len(seq)

        yield seq
    
def convert_to_fastq_if_needed(filepath, output_dir, tmpfilename='temp.fa'):
    '''
    Convert a FASTA file provided as input to a FASTQ if needed
    '''
    # if it is a FASTA file the first character is >
    with open(filepath, 'r') as f:
        is_fasta = f.read(1) == '>'
    
    if not is_fasta:
        return filepath
    else:
        SeqIO.write(
            fasta_to_fastq_record(SeqIO.parse(filepath, 'fasta')),
            os.path.join(output_dir, tmpfilename), 'fastq'
        )
        return os.path.join(output_dir, tmpfilename)

def load_and_validate_input(input_yaml, output_dir):
    '''
    Check for errors in the small RNA pipeline input YAML
    '''
    # TODO: check file types and existance
    with open(input_yaml) as f:
        input_yaml = yaml.load(f.read(), Loader=yaml.Loader)

    smallRNA_fastq_present = 'smallRNA_fastq' in input_yaml
    size_sorted_fastqs_present = 'size_sorted_fastqs' in input_yaml

    # smallRNA_fastq and size_sorted_fastqs are mutaually exclusive in the input file
    if (smallRNA_fastq_present and size_sorted_fastqs_present) or not (smallRNA_fastq_present or size_sorted_fastqs_present):
        return False, 'Exactly one of smallRNA_fastq or size_sorted_fastqs must be specified in the config file', input_yaml

    if smallRNA_fastq_present:
        input_yaml['smallRNA_fastq'] = convert_to_fastq_if_needed(
            unzip_if_needed(input_yaml['smallRNA_fastq'], output_dir, tmpfilename='smallRNA.fastq'), 
            output_dir, tmpfilename='smallRNA.fasta.fastq'
        )

    # cds and unspliced_transcriptome set to None if not present
    if 'cds' not in input_yaml:
        input_yaml['cds'] = None
    else:
        input_yaml['cds'] = unzip_if_needed(input_yaml['cds'], output_dir, tmpfilename='cds.fasta')

    if 'unspliced_transcriptome' not in input_yaml:
        input_yaml['unspliced_transcriptome'] = None
    else:
        input_yaml['unspliced_transcriptome'] = unzip_if_needed(input_yaml['unspliced_transcriptome'], output_dir, tmpfilename='unspliced.fasta')

    # program paths are all optional, default to in path
    if 'program_paths' not in input_yaml:
        input_yaml['program_paths'] = {}

    if 'unitas' not in input_yaml['program_paths']:
        input_yaml['program_paths']['unitas'] = 'unitas.pl'

    if 'bowtie2-build' not in input_yaml['program_paths']:
        input_yaml['program_paths']['bowtie2-build'] = 'bowtie2-build'

    if 'bowtie2' not in input_yaml['program_paths']:
        input_yaml['program_paths']['bowtie2'] = 'bowtie2'

    if 'cutadapt' not in input_yaml['program_paths']:
        input_yaml['program_paths']['cutadapt'] = 'cutadapt'

    if 'samtools' not in input_yaml['program_paths']:
        input_yaml['program_paths']['samtools'] = 'samtools'

    if 'bedtools' not in input_yaml['program_paths']:
        input_yaml['program_paths']['bedtools'] = 'bedtools'

    # trim
    if 'trim' in input_yaml:
        # if trim is specified, require smallRNA_fastq
        if not smallRNA_fastq_present:
            return False, 'For trim and sort steps smallRNA_fastq must be specified, not size_sorted_fastqs', input_yaml
        else:
            input_yaml['trim']['input'] = input_yaml['smallRNA_fastq']

        # kit is mually exclusive to 5_prime and 3_prime
        if ('5_prime' in input_yaml['trim'] or '3_prime' in input_yaml['trim']) and 'kit' in input_yaml['trim']:
            return False, 'If a kit is specified in the trim config, 5_prime and 3_prime sequences should not be', input_yaml
        
        # kit must be one of qiagen, qiagen-rev, neb or neb-rev
        if 'kit' in input_yaml['trim'] and input_yaml['trim']['kit'] not in KNOWN_ADAPTERS.keys():
            return False, f'The kit specified in trim should be one of: {", ".join(KNOWN_ADAPTERS.keys())}', input_yaml
        elif 'kit' in input_yaml['trim']:
            input_yaml['trim']['5_prime'] = KNOWN_ADAPTERS[input_yaml['trim']['kit']]['5']
            input_yaml['trim']['3_prime'] = KNOWN_ADAPTERS[input_yaml['trim']['kit']]['3']

            del input_yaml['trim']['kit']
        
        # min_quality set to 20 if it doesn't exist
        if 'min_quality' not in input_yaml['trim']:
            input_yaml['trim']['min_quality'] = 20
     
    # sort
    if 'sort' in input_yaml:
        if not smallRNA_fastq_present:
            return False, 'For trim and sort steps smallRNA_fastq must be specified, not size_sorted_fastqs', input_yaml
        elif 'trim' in input_yaml:
            input_yaml['sort']['input'] = os.path.join(output_dir, 'trimmed_reads.fq')
        else:
            input_yaml['sort']['input'] = input_yaml['smallRNA_fastq']

        if 'genome' in input_yaml['sort']:
            input_yaml['sort']['genome'] = unzip_if_needed(input_yaml['sort']['genome'], output_dir, tmpfilename='genome.fasta')

        input_yaml['sort']['cds'] = input_yaml['cds']

        # align_to_cds defaults to false, if true cds must be present
        if 'align_to_cds' not in input_yaml['sort']:
            input_yaml['sort']['align_to_cds'] = False
        elif input_yaml['sort']['align_to_cds'] and input_yaml['sort']['cds'] is None:
            return False, 'If align_to_cds is specified in the sort section, a cds file must be specified as well', input_yaml
        
        # min_length defaults to 15
        if 'min_length' not in input_yaml['sort']:
            input_yaml['sort']['min_length'] = 15
        
        # max_length defaults to 30
        if 'max_length' not in input_yaml['sort']:
            input_yaml['sort']['max_length'] = 30

        # mismatches defaults to 0
        if 'mismatches' not in input_yaml['sort']:
            input_yaml['sort']['mismatches'] = 0
        
    # unitas
    if 'unitas' in input_yaml:
        # if sort is not run size_sorted_fastqs must be specified 
        if 'sort' not in input_yaml and 'size_sorted_fastqs' not in input_yaml:
            return False, 'To run unitas, either the sort step must be run first or size_sorted_fastqs provided instead of smallRNA_fastq', input_yaml
        elif 'size_sorted_fastqs' in input_yaml:
            input_yaml['unitas']['inputs'] = input_yaml['size_sorted_fastqs']
        elif 'sort' in input_yaml:
            input_yaml['unitas']['inputs'] = os.path.join(output_dir, 'binned_length_rna')

        # at least one of refseq or species must be specified 
        if 'refseq' not in input_yaml['unitas'] and ('species' not in input_yaml['unitas'] or input_yaml['unitas']['species'] == 'x'):
            return False, 'At least one refseq or non-x species must be provided in the unitas config, to reuse cds and unsplices transcriptiome, specify gene', input_yaml
        
        if 'species' not in input_yaml['unitas']:
            input_yaml['unitas']['species'] = 'x'

        if 'refseq' in input_yaml['unitas'] and 'gene' in input_yaml['unitas']['refseq']:
            gene_index = input_yaml['unitas']['refseq'].index('gene')
            del input_yaml['unitas']['refseq'][gene_index]

            if input_yaml['cds'] is None and input_yaml['unspliced_transcriptome'] is None:
                return False, 'If gene is specified in unitas or targetid refseq or unspliced_transcriptome must also be specified', input_yaml
            
            if input_yaml['cds'] is not None:
                input_yaml['unitas']['refseq'].append({'gene': input_yaml['cds']})
            
            if input_yaml['unspliced_transcriptome'] is not None:
                input_yaml['unitas']['refseq'].append({'gene': input_yaml['unspliced_transcriptome']})
        
    # targetid
    if 'targetid' in input_yaml:
        if 'sort' in input_yaml:
            input_yaml['targetid']['inputs'] = os.path.join(output_dir, 'binned_length_rna')
        elif 'size_sorted_fastqs' in input_yaml:
            input_yaml['targetid']['inputs'] = input_yaml['size_sorted_fastqs']
        else:
            input_yaml['targetid']['inputs'] = input_yaml['smallRNA_fastq']

        # at least one target file specified
        if 'target_files' not in input_yaml['targetid'] or len(input_yaml['targetid']['target_files']) < 1:
            return False, 'To use the targetid command, you must specify at least one target file, to reuse cds and unsplices transcriptiome, specify gene', input_yaml

        if 'gene' in input_yaml['targetid']['target_files']:
            gene_index = input_yaml['targetid']['target_files'].index('gene')
            del input_yaml['targetid']['target_files'][gene_index]

            if input_yaml['cds'] is None and input_yaml['unspliced_transcriptome'] is None:
                return False, 'If gene is specified in unitas or targetid refseq or unspliced_transcriptome must also be specified', input_yaml
            
            if input_yaml['cds'] is not None:
                input_yaml['targetid']['target_files'].append(input_yaml['cds'])
            
            if input_yaml['unspliced_transcriptome'] is not None:
                input_yaml['targetid']['target_files'].append(input_yaml['unspliced_transcriptome'])

        # min_seq_length defaults to 5
        if 'min_seq_length' not in input_yaml['targetid']:
            input_yaml['targetid']['min_seq_length'] = 5

        # mismatches defaults to 0
        if 'mismatches' not in input_yaml['targetid']:
            input_yaml['targetid']['mismatches'] = 0

    return True, None, input_yaml