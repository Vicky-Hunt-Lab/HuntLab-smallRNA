import os
import os.path
import glob
import csv
import re

from subprocess import run
from collections import defaultdict

import numpy as np

from matplotlib import pyplot as plt

from .config import do_log, get_config_key

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
    file_dict_int = defaultdict(lambda: [])
    file_dict_str = defaultdict(lambda: [])
    new_file = [[]]

    for filename in glob.glob(os.path.join(get_config_key('general', 'output_directory'), 'unitas', '*/unitas.annotation_summary.txt')):
        with open(filename) as f:
            reader = csv.reader(f, delimiter='\t')

            length = re.findall(r'length(\d+)', filename)
            if len(length) < 1:
                file_table = [['File Name : ', filename]]
            else:
                file_table = [['RNA Length ' + length[0], '']]

            for line in reader:
                file_table.append(line)

            if len(length) > 0:
                file_dict_int[int(length[0])] = file_table
            else:
                file_dict_str[str(filename)] = file_table

    for key in sorted(file_dict_int.keys()):
        file_table = file_dict_int[key]

        while len(new_file) < len(file_table):
            new_file.append(['' for i in range(len(new_file[0]))])

        for i, line in enumerate(file_table):
            new_file[i] = new_file[i] + line + ['']

    for key in sorted(file_dict_str.keys()):
        file_table = file_dict_str[key]

        while len(new_file) < len(file_table):
            new_file.append(['' for i in range(len(new_file[0]))])

        for i, line in enumerate(file_table):
            new_file[i] = new_file[i] + line + ['']

    unitas_table = os.path.join(get_config_key('general', 'output_directory'), 'unitas_summery.csv')
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

        headers = next(reader)
        values = defaultdict(lambda: {})

        data_rows = list(reader)

        fig, ax = plt.subplots(figsize=(11.69, 8.27))

        labels = set()

        for i, header in enumerate(headers):
            if header != '':
                rna_length = re.findall(r'\d+', header)

                if len(rna_length) > 0:
                    rna_length = rna_length[0]

                    for row in data_rows:
                        try:
                            if len(row[i]) > 0 and not row[i][0].isspace():
                                values[rna_length][row[i]] = float(row[i + 1])
                        except ValueError:
                            pass
                        except IndexError:
                            break

        for length in values:
            labels = labels | set(values[length].keys())

        with open(os.path.join(get_config_key('general', 'output_directory'), 'unitas_graph_data.csv'), 'w') as f:
            writer = csv.writer(f)

            rna_lengths = list(values.keys())
            writer.writerow(['RNA Length'] + rna_lengths)
        
            offsets = None
            for l in labels:
                counts = []

                for key in values.keys():
                    try:
                        counts.append(values[key][l])
                    except KeyError:
                        counts.append(0)

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