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
import glob
import csv
import re

from subprocess import run
from collections import defaultdict

import numpy as np

from matplotlib import pyplot as plt
from Bio import SeqIO

from .config import do_log, get_config_key

def copy_and_label_file(ref_file, filename, label):
    '''
    Labeling function for the CDS and Unspliced Transcriptome files
    '''

    def rename_iter(path):
        for seq in SeqIO.parse(path, 'fasta'):
            seq.id = f'{label}|{seq.id}'

            yield seq

    
    new_path = os.path.join(get_config_key('general', 'output_directory'), f'{filename}.fasta')
    SeqIO.write(rename_iter(ref_file), new_path, 'fasta')

    return new_path

def run_unitas_annotation(small_rna, species_name, ref_files, quiet=0, unitas_output='.'):
    '''
    Run Unitas on the small RNAs, given a file and species name
    '''
    CWD = os.getcwd()
    SMALL_RNA_PATH = os.path.join(CWD, small_rna)

    os.chdir(unitas_output)

    unitas_command = [
        get_config_key('cli-tools', 'unitas', 'path_to_unitas'),
        '-input', SMALL_RNA_PATH,
        '-species', species_name
    ]

    for refs in ref_files:
        unitas_command.append('-refseq')
        unitas_command.append(os.path.join(CWD, refs))

    unitas_command = unitas_command + get_config_key('cli-tools', 'unitas', 'unitas_params')

    if get_config_key('cli-tools', 'unitas', 'unitas_pass_threads'):
        threads = get_config_key('general', 'threads')

        unitas_command.append('-threads')
        unitas_command.append(str(threads))

    do_log(quiet, '==> Running first unitas pass')
    result = run(unitas_command, capture_output=(quiet != 0))

    result.check_returncode()

    os.chdir(CWD)

def merge_summary():
    '''
    Merge together the summery files into one CSV
    '''
    new_file = []

    for filename in glob.glob(os.path.join(get_config_key('general', 'output_directory'), 'unitas', '*/unitas.annotation_summary.txt')):
        with open(filename) as f:
            reader = csv.reader(f, delimiter='\t')

            length = re.findall(r'length(\d+)', filename)
            if len(length) < 1:
                new_file.append(['File Name : ', filename])
            else:
                new_file.append(['RNA Length ' + length[0], ''])

            names = []
            values = []
            for line in reader:
                names.append(line[0])
                values.append(line[1])

            new_file.append(names)
            new_file.append(values)

    unitas_table = os.path.join(get_config_key('general', 'output_directory'), 'unitas_summary.csv')
    with open(unitas_table, 'w') as f:
        writer = csv.writer(f)

        writer.writerows(new_file)

    return unitas_table

def graph_unitas_classification_type(path_to_table):
    '''
    Craete a graph of small RNA by length against unitas type
    '''
    with open(path_to_table) as f:
        reader = csv.reader(f)

        fig, ax = plt.subplots(figsize=(11.69, 8.27))

        labels = set()
        values = defaultdict(lambda: defaultdict(lambda: 0))

        # Use skip_group to skip any group not named by length
        skip_group = False
        for i, line in enumerate(reader):
            if i % 3 == 0:
                if not line[0].startswith('RNA Length '):
                    skip_group = True
                else:
                    line_name = re.findall('\d+', line[0])[0]
            elif i % 3 == 1:
                if skip_group:
                    continue
                
                line_labels = set()
                for l in line:
                    if not l[0].isspace():
                        line_labels.add(l)

                all_labels = line
                        
                labels = labels | line_labels
            elif i % 3 == 2:
                if skip_group:
                    skip_group = False
                    continue

                for j, name in enumerate(all_labels):
                    if name in line_labels:
                        values[line_name][name] = float(line[j])

        with open(os.path.join(get_config_key('general', 'output_directory'), 'unitas_graph_data.csv'), 'w') as f:
            writer = csv.writer(f)

            rna_lengths = list(values.keys())
            writer.writerow(['RNA Length'] + rna_lengths)
        
            offsets = None
            for l in labels:
                counts = []

                for key in rna_lengths:
                    counts.append(values[key][l])

                ax.bar(rna_lengths, counts, label=l, bottom=offsets)

                writer.writerow([l] + counts)

            if offsets is None:
                offsets = np.array(counts)
            else:
                offsets = offsets + np.array(counts)

    ax.set_xticklabels(rna_lengths, fontsize=7)
    ax.set_xlabel('Length of small RNA')
    ax.set_ylabel('')

    plt.legend(title='Unitas Category')

    plt.savefig(os.path.join(get_config_key('general', 'output_directory'), 'unitasGraph.png'))

if __name__ == '__main__':
    graph_unitas_classification_type('output/unitas_summery.csv')