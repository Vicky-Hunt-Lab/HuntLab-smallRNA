import toml

import os

from multiprocessing import cpu_count

CONFIG = {}

# TODO: thread passing to programs

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
            'path_to_bbmap': 'bbmap.sh',
            'bbmap_align_params': [],
            'bbmap_index_params': []
        },
        'unitas': {
            'path_to_unitas': 'unitas',
            'unitas_params': []
        }
    }
}

def load_config(path_to_config):
    '''
    Load the config file into the global CONFIG varible and validate the result
    '''
    global CONFIG
    CONFIG = toml.load(path_to_config)

    result = validate_toml(CONFIG, DEFALT_CONFIG)

    if result is not None:
        raise Exception(result)

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