from subprocess import run

from config import get_config_key

def run_unitas_annotation(small_rna, species_name):
    '''
    Run Unitas on the small RNAs, given a file and species name
    '''
    unitas_command = [
        get_config_key('cli-tools', 'unitas', 'path_to_unitas'),
        '-input', small_rna,
        '-species', species_name
    ]

    unitas_command = unitas_command + get_config_key('cli-tools', 'unitas', 'unitas_params')

    run(unitas_command)