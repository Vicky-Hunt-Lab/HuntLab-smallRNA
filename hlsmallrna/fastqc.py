import os.path
import zipfile
import glob
import csv

from subprocess import run
from io import StringIO
from math import inf

from Bio import SeqIO

from .config import get_config_key, mkdir_if_not_exists, do_log

def run_fastqc(path_to_fastq, quiet=0):
    '''
    Call fastqc and trim off any bases with a lower quatile below 20
    '''
    outdir = os.path.join(get_config_key('general', 'output_directory'), 'fastqc')
    mkdir_if_not_exists(outdir)

    command = [
        get_config_key('cli-tools', 'fastqc', 'path_to_fastqc'), 
        '--outdir', outdir, 
        path_to_fastq
    ]

    command = command + get_config_key('cli-tools', 'fastqc', 'fastqc_params')

    if get_config_key('cli-tools', 'fastqc', 'fastqc_pass_threads'):
        threads = get_config_key('general', 'threads')

        command.append('-t')
        command.append(str(threads))

    do_log(quiet, '====> Running FastQC')
    result = run(command, capture_output=(quiet != 0))
    result.check_returncode()

    return result.returncode == 0

def cut_rna_below_cutoff(fastq_path, cutoff, quiet=0):
    '''
    Read the fastqc result and automatically trim off any base that is below the threshold
    '''

    do_log(quiet, '====> Extracting results from FastQC Output')

    # Unzip file produced by fastqc
    path_to_zip = os.path.join(get_config_key('general', 'output_directory'), 'fastqc')
    zip_file = glob.glob(os.path.join(path_to_zip, '*.zip'))[0]
    
    with zipfile.ZipFile(zip_file, 'r') as f:
        f.extractall(path=path_to_zip)

    # parse in quality score data, cut out part of file
    with open('.' + os.path.join(zip_file.strip('.zip'), 'fastqc_data.txt')) as f:
        data = f.readlines()

    start_line = None
    end_line = None
    for i, line in enumerate(data):
        if line.startswith('>>Per base sequence quality'):
            start_line = i + 2

        if start_line is not None and line.startswith('>>END_MODULE'):
            end_line = i
            break

    # use csv reader for simplicity
    cutoff_base = inf

    reader = reversed(list(csv.reader(StringIO(''.join(data[start_line:end_line])), delimiter='\t')))
    for line in reader:
        if float(line[3]) < cutoff:
            if cutoff_base == int(line[0]) + 1:
                cutoff_base = int(line[0])
        elif cutoff_base == inf:
            if float(line[3]) < cutoff:
                cutoff_base = int(line[0])
            else:
                cutoff_base = None

    # cut sequences from 3' end
    cut_func = lambda x: x[:cutoff_base]

    sequences = SeqIO.parse(fastq_path, 'fastq')
    cut_seqs = map(cut_func, sequences)
    SeqIO.write(cut_seqs, os.path.join(get_config_key('general', 'output_directory'), 'cut_sequences.fastq'), 'fastq')