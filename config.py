import toml

from multiprocessing import cpu_count

try:
    CONFIG = toml.load('config.toml')
except:
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
            # TODO: Put an actual command here
            'path_to_trim': 'trim',
            'trim_type': '',
            'trim_params': []
        },
        'fastqc': {
            'fastqc_pass_threads': True,
            'path_to_fastqc': 'fastqc',
            'fastqc_params': []
        },
        'bowtie2': {
            'path_to_bowtie2': 'bowtie2',
            'path_to_bowtie2_build': 'bowtie2-build',
            'bowtie2_params': [],
            'bowtie2_build_params': []
        },
        'samtools': {
            'path_to_samtools': 'samtools',
            'samtools_view_params': [],
            'samtools_sort_params': []
        },
        'bedtools': {
            'path_to_bedtools': 'bedtools',
            'bedtools_bamToFastq_params': []
        },
        'unitas': {
            'path_to_unitas': 'unitas',
            'unitas_params': []
        }
    }
}

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