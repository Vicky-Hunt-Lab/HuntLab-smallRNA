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
import glob
import csv
import re

from math import inf
from subprocess import run
from collections import defaultdict

import numpy as np

from matplotlib import pyplot as plt
from Bio import SeqIO

def find_min_max_lengths(unitas_in):
    '''
    Calculate the minimum and maximum lengths for the unitas report based on fastq filenames
    '''
    minimum = inf
    maximum = -inf

    filenames = glob.glob(os.path.join(unitas_in, '*.fq')) + glob.glob(os.path.join(unitas_in, '*.fastq'))
    for name in filenames:
        length = int(re.findall(r'length(\d+)', name)[0])

        if length > maximum:
            maximum = length
        
        if length < minimum:
            minimum = length

    if maximum == -inf:
        maximum = 30
    if minimum == inf:
        minimum = 18

    return minimum, maximum

def run_unitas_annotation(species_name, ref_files, program_paths, output, verbose, small_rna):
    '''
    Run Unitas on the small RNAs, given a file and species name
    '''
    CWD = os.getcwd()
    SMALL_RNA_PATH = os.path.join(CWD, small_rna)

    os.chdir(output)

    unitas_command = [
        program_paths['unitas'],
        '-input', SMALL_RNA_PATH,
        '-species', species_name
    ]

    if ref_files is not None:
        for refs in ref_files:
            unitas_command.append('-refseq')
            unitas_command.append(os.path.join(CWD, refs))

    print(f'====> Running unitas on {os.path.basename(small_rna)}...')
    result = run(unitas_command, capture_output=not verbose)

    result.check_returncode()

    os.chdir(CWD)

def merge_summary(output):
    '''
    Merge together the summery files into one CSV
    '''
    numbered_names = defaultdict(lambda: 'Unnamed')
    numbered_results = defaultdict(lambda: defaultdict(lambda: 0))
    named_results = defaultdict(lambda: defaultdict(lambda: 0))
    unitas_order = defaultdict(lambda: [])

    for filename in glob.glob(os.path.join(output, 'unitas', '*/unitas.annotation_summary.txt')):
        length_regex_result = re.findall(r'length(\d+)', filename)
        if len(length_regex_result) < 1:
            namestr = f'File Name : {filename}'
        else:
            length = int(length_regex_result[0])
            namestr = f'RNA Length {length}'
            numbered_names[length] = namestr

            namestr = length

        with open(filename) as f:
            reader = csv.reader(f, delimiter='\t')

            for line in reader:
                if type(namestr) == int:
                    numbered_results[namestr][line[0]] = float(line[1])
                else:
                    named_results[namestr][line[0]] = float(line[1])

                unitas_order[namestr].append(line[0])

    combined_unitas_file = []
    for result in sorted(numbered_results.keys()):
        names = []
        values = []

        for unitas_cat in unitas_order[result]:
            names.append(unitas_cat)
            values.append(numbered_results[result][unitas_cat])

        combined_unitas_file.append([numbered_names[result]])
        combined_unitas_file.append(names)
        combined_unitas_file.append(values)
        combined_unitas_file.append([])

    for result in sorted(named_results.keys()):
        names = []
        values = []

        for unitas_cat in unitas_order[result]:
            names.append(unitas_cat)
            values.append(numbered_results[result][unitas_cat])

        combined_unitas_file.append([result])
        combined_unitas_file.append(names)
        combined_unitas_file.append(values)
        combined_unitas_file.append([])

    unitas_table = os.path.join(output, 'unitas_summary.csv')
    with open(unitas_table, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(combined_unitas_file)
    
    return numbered_results

def graph_unitas_classification_type(numbered_results, min_length, max_length, output):
    '''
    Craete a graph of small RNA by length against unitas type
    '''
    fig, ax = plt.subplots(figsize=(11.69, 8.27))
    final_cats = ['gene', '5_UTR', '3_UTR', 'miRNA', 'tRNA', 'rRNA', 'TE', 'no annotation', 'low_complexity']
    rna_lengths = list(range(min_length, max_length + 1))

    # figure out if there are any other toplevel categories we need to include on the graph
    included_cats = set(final_cats)
    for length in numbered_results:
        for category in numbered_results[length]:
            if not category.startswith(' '):
                included_cats.add(category)

    for label in included_cats:
        if label not in final_cats:
            final_cats.append(label)

    # Plot each caegory on the graph and add to the plot data table
    base_offset = np.array([0] * len(rna_lengths))
    plot_data = [['RNA Length'] + rna_lengths]
    for category in final_cats:
        category_bars = []
        for length in rna_lengths:
            if category in numbered_results[length]:
                category_bars.append(numbered_results[length][category])
            else:
                category_bars.append(0)

        ax.bar(rna_lengths, category_bars, bottom=base_offset, label=category, align='edge', width=0.8)
        plot_data.append([category] + category_bars)
        base_offset = base_offset + np.array(category_bars)

    with open(os.path.join(output, 'unitas_graph_data.csv'), 'w') as f:
        writer = csv.writer(f)
        writer.writerows(plot_data)

    ax.set_xticks(rna_lengths)
    ax.set_xticklabels(rna_lengths, fontsize=7)
    ax.set_xlabel('Length of small RNA')
    ax.set_ylabel('Unitas RPM')

    plt.legend(title='Unitas Category')

    plt.savefig(os.path.join(output, 'unitasGraph.png'))