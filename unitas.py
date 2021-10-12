import os
import os.path
import glob
import csv
import re

from subprocess import run
from collections import defaultdict

from config import get_config_key

def run_unitas_annotation(small_rna, species_name, ref_files, quiet=False, unitas_output='.'):
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

    print('==> Running first unitas pass')
    result = run(unitas_command, capture_output=quiet)

    result.check_returncode()

    os.chdir(CWD)

    # TODO: second pass of unitas - how best to do that?

def merge_summary():
    '''
    Merge together the summery files into one CSV
    '''
    file_dict = defaultdict(lambda: [])
    new_file = [[]]

    for filename in glob.glob(os.path.join(get_config_key('general', 'output_directory'), 'unitas', '*/unitas.annotation_summary.txt')):
        with open(filename) as f:
            reader = csv.reader(f, delimiter='\t')

            length = re.findall(r'length(\d+)', filename)
            if len(length) < 1:
                file_table = [['File Name : ', filename]]
            else:
                file_table = [['RNA Length ' + length[0], '']]

            te_count = 0
            for line in reader:
                file_table.append(line)
                if len(file_table) >= 28 and 'no annotation' not in line[0]:
                    te_count += float(line[1])
                    file_table[-1][0] = '   ' + file_table[-1][0]

            file_table = file_table[:27] + [['TEs', te_count]] + file_table[27:]

            if len(length) > 0:
                file_dict[int(length[0])] = file_table
            else:
                file_dict[filename] = file_table

    for key in sorted(file_dict.keys()):
        file_table = file_dict[key]

        while len(new_file) < len(file_table):
            new_file.append(['' for i in range(len(new_file[0]))])

        for i, line in enumerate(file_table):
            new_file[i] = new_file[i] + line + ['']

    with open(os.path.join(get_config_key('general', 'output_directory'), 'unitas_summery.csv'), 'w') as f:
        writer = csv.writer(f)

        writer.writerows(new_file)
