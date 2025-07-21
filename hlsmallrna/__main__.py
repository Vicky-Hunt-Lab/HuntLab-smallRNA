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
import glob
import shutil

from functools import partial
from multiprocessing import Pool
from argparse import ArgumentParser

from .validate_config import load_and_validate_input
from .trim import run_trim
from .genome_align import align_to_genome, bin_rna_size, graph_length
from .label_for_unitas import label_file_for_unitas
from .unitas import find_min_max_lengths, run_unitas_annotation, merge_summary, graph_unitas_classification_type
from .targetid import merge_and_revcomp_input_files, find_targets, build_summary_files, do_term_enrichment

def run_trim_reads(trim_config, program_paths, output, threads, verbose, keep_intermediate):
    '''
    Call run_trim with the correct values based on the config
    '''
    run_trim(
        trim_config['input'], output, five_prime_adapter=trim_config['5_prime'], three_prime_adapter=trim_config['3_prime'], quality_cutoff=trim_config['min_quality'], 
        threads=threads, cutadapt_path=program_paths['cutadapt'], keep_intermediate=keep_intermediate, verbose=verbose
    )

def run_sort_reads(sort_config, program_paths, output, threads, verbose, keep_intermediate):
    '''
    Call align_to_genome with the correct values based on the config
    '''
    if 'genome' in sort_config:
        sort_input = align_to_genome(
            sort_config['genome'], sort_config['input'], sort_config['cds'], program_paths, output,
            verbose=verbose, keep_intermediate=keep_intermediate, threads=threads, mismatches=sort_config['mismatches']
        )
    else:
        print('====> No genome specified, skipping align...')
        sort_input = sort_config['input']

    table_file = bin_rna_size(sort_input, sort_config['min_length'], sort_config['max_length'], output, verbose=verbose)
    graph_length(table_file, output)

def run_unitas_on_lenths(unitas_config, min_length, max_length, program_paths, output, threads, verbose, keep_intermediate):
    '''
    Call label_for_unitas and unitas commands based on the config file
    '''
    REFSEQ_LABELED_DIR = os.path.join(output, 'unitas_refseq_labelled')
    UNITAS_OUTPUT_DIR = os.path.join(output, 'unitas')

    if not os.path.exists(REFSEQ_LABELED_DIR):
        os.makedirs(REFSEQ_LABELED_DIR)

    if not os.path.exists(UNITAS_OUTPUT_DIR):
        os.makedirs(UNITAS_OUTPUT_DIR)

    # First sort out any refseq files specified
    refseq_paths = []

    if 'refseq' in unitas_config:
        for refseq in unitas_config['refseq']:
            if type(refseq) == str:
                refseq_paths.append(refseq)
            elif type(refseq) == dict:
                for k in refseq.keys():
                    labelled_file = os.path.join(REFSEQ_LABELED_DIR, f'{os.path.basename(refseq[k])}.labelled.fasta')
                    label_file_for_unitas(refseq[k], labelled_file, k)

                    refseq_paths.append(labelled_file)

    # build and run unitas functions in parrellel
    small_rna_files = glob.glob(os.path.join(unitas_config['inputs'], '*.fq')) + glob.glob(os.path.join(unitas_config['inputs'], '*.fastq'))

    unitas_mappable_function = partial(run_unitas_annotation, unitas_config['species'], refseq_paths, program_paths, UNITAS_OUTPUT_DIR, verbose)
    with Pool(threads) as p:
        p.map(unitas_mappable_function, small_rna_files)

    unitas_length_results = merge_summary(output)
    graph_unitas_classification_type(unitas_length_results, min_length, max_length, output)

    if not keep_intermediate:
        shutil.rmtree(os.path.join(output, 'unitas_refseq_labelled'))

def run_targetid(targetid_config, program_paths, output, threads, verbose, keep_intermediate):
    '''
    Call the targetid functions based on the settings in the config file
    '''
    if os.path.isdir(targetid_config['inputs']):
        small_rna_input = glob.glob(os.path.join(targetid_config['inputs'], '*.fq')) + glob.glob(os.path.join(targetid_config['inputs'], '*.fastq'))
    else:
        small_rna_input = [targetid_config['inputs']]

    revcomp_file = merge_and_revcomp_input_files(small_rna_input, output)
    target_sam_files = find_targets(
        revcomp_file, targetid_config['target_files'], output, program_paths, threads=threads, 
        min_seq_length=targetid_config['min_seq_length'], mismatches_allowed=targetid_config['mismatches'],
        verbose=verbose
    )
    targets = build_summary_files(target_sam_files, output, program_paths, verbose=verbose)

    if 'enrich' in targetid_config:
        do_term_enrichment(
            targetid_config['enrich']['included_files'], targets, targetid_config['enrich']['eggnog_data'], output, program_paths,
            threads=threads, verbose=verbose
        )

    if not keep_intermediate:
        shutil.rmtree(os.path.join(output, 'bowtie_indexes'))
        shutil.rmtree(os.path.join(output, 'target_alignments'))
        os.remove(os.path.join(output, 'revcomp_rna.fastq'))

def main():
    parser = ArgumentParser(description='Pipeline to process small RNAs - version 2')

    parser.add_argument('config_file', help='Path to YAML config file')

    parser.add_argument('-o', '--output', help='Directory to write output to (will be created if it doesn\'t exist, default: hlsmallrna_output)', default='hlsmallrna_output')
    parser.add_argument('-t', '--threads', help='Number of threads to run the pipeline on (default: 4)', type=int, default=4)
    parser.add_argument('-k', '--keep-files', help='If set, keep the intermediate files, to help debug the application', action='store_true')
    parser.add_argument('-v', '--verbose', help='If set, print the output for the intermediate commands', action='store_true')

    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    shutil.copy2(args.config_file, os.path.join(args.output, 'config.yml'))

    print('==> Validating config and converting files...')
    success, err_msg, config = load_and_validate_input(args.config_file, args.output)

    if not success:
        print(f'Error: {err_msg}')
        exit(1)
    
    if 'trim' in config:
        print('==> Trimming Reads...')
        run_trim_reads(config['trim'], config['program_paths'], args.output, args.threads, args.verbose, args.keep_files)

    if 'sort' in config:
        print('==> Aligning to genome and sorting reads by length...')
        run_sort_reads(config['sort'], config['program_paths'], args.output, args.threads, args.verbose, args.keep_files)
        
    if 'unitas' in config:
        if 'sort' in config:
            min_length = config['sort']['min_length']
            max_length = config['sort']['max_length']
        else:
            min_length, max_length = find_min_max_lengths(config['unitas']['inputs'])

        print('==> Running Unitas in parellel on each sequence length...')
        run_unitas_on_lenths(config['unitas'], min_length, max_length, config['program_paths'], args.output, args.threads, args.verbose, args.keep_files)

    if 'targetid' in config:
        print('==> Finding targets of small RNA...')
        run_targetid(config['targetid'], config['program_paths'], args.output, args.threads, args.verbose, args.keep_files)