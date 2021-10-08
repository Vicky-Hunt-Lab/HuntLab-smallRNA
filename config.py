from posixpath import pardir
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
    global CONFIG
    CONFIG = toml.load(path_to_config)

def get_config_key(*keys):
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
    if not os.path.exists(path):
        os.mkdir(path)