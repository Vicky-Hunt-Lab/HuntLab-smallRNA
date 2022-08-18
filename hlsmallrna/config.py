import toml

import os

from multiprocessing import cpu_count

# global object to load config into, it works, but isn't nice
# or good practice
CONFIG = {}

DEFALT_CONFIG = {
    'general': {
        'threads': cpu_count(),
        'output_directory': './output'
    },
    'cli-tools': {
        'trim': {
            'trim_pass_threads': True,
            'path_to_trim': 'cutadapt',
            'trim_type': '',
            'trim_params': []
        },
        'fastqc': {
            'fastqc_pass_threads': True,
            'path_to_fastqc': 'fastqc',
            'fastqc_params': []
        },
        'bbmap': {
            'bbmap_pass_threads': True,
            'path_to_bbmap': 'bbmap.sh',
            'bbmap_align_params': [],
            'bbmap_index_params': []
        },
        'unitas': {
            'unitas_pass_threads': True,
            'path_to_unitas': 'unitas.pl',
            'unitas_params': []
        },
        'bowtie2': {
            'bowtie2_pass_threads': True,
            'path_to_bowtie2': 'bowtie2',
            'bowtie2_params': [],
            'path_to_bowtie2_build': 'bowtie2-build',
            'bowtie2_build_params': []
        },
        'samtools': {
            'path_to_samtools': 'samtools',
            'samtools_view_params': ['-h', '-F', '256', '-F', '4'],
            'samtools_fastq_params': []
        },
        'bedtools': {
            'path_to_bedtools': 'bedtools',
            'bedtools_bamtofastq_params': []
        }
    }
}

def load_config(path_to_config, quiet=0):
    '''
    Load the config file into the global CONFIG varible and validate the result
    '''
    global CONFIG

    try:
        CONFIG = toml.load(path_to_config)

        result = validate_toml(CONFIG, DEFALT_CONFIG)

        if result is not None:
            raise Exception(result)
    except FileNotFoundError:
        do_log(quiet, f'No {path_to_config} found, running with defaults...')

def validate_toml(in_dict, compare):
    '''
    Check the type in a TOML dict (in_dict) match the types of another dict,
    validates the types in a loaded config file
    '''
    for key in compare:
        if key in in_dict:
            if type(in_dict[key]) != type(compare[key]):
                print('Error in config file!')

                return key + ' should be ' + str(type(compare[key])) + ' but found ' + str(type(in_dict[key]))

            if type(in_dict[key]) == dict:
                result = validate_toml(in_dict[key], compare[key])

                if result is not None:
                    return result

def get_config_key(*keys):
    '''
    Retrieve a config key from the config file if it exists or the
    default config if it doesn't
    '''
    try:
        last_item = CONFIG
        
        for key in keys:
            last_item = last_item[key]
    except KeyError:
        last_item = DEFALT_CONFIG

        for key in keys:
            last_item = last_item[key]

    return last_item

def mkdir_if_not_exists(path):
    '''
    Shorthand to make directories only if they don't exist already
    without crashing
    '''
    if not os.path.exists(path):
        os.mkdir(path)

def do_log(quiet, *print_args, **print_kwds):
    '''
    Function that prints if quiet isn't set
    '''
    if quiet < 2:
        print(*print_args, **print_kwds)