import os

from subprocess import run

from .config import get_config_key, mkdir_if_not_exists, do_log

def samestrand_overlap(genome_file, rna_file_1, rna_file_2, longest_rna, quiet=0):
    '''
    Run the samestand overlap script
    '''

    do_log(quiet, '==> Running same strand overlap script')

    CWD = os.getcwd()
    OUTPUT_DIR = os.path.join(CWD, get_config_key('general', 'output_directory'), 'samestrand_overlap')

    mkdir_if_not_exists(OUTPUT_DIR)
    os.chdir(OUTPUT_DIR)

    command = [
        'overlap_ss.sh',
        os.path.join(CWD, genome_file),
        os.path.join(CWD, rna_file_1),
        os.path.join(CWD, rna_file_2),
        str(longest_rna)
    ]

    run(command, capture_output=(quiet != 0))

    os.chdir(CWD)

    do_log(quiet, '==> Finished same strand overlap script')