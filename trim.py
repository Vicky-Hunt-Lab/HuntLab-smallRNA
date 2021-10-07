import os.path

from subprocess import run

from config import get_config_key

def run_trim(file_to_trim, start_adapter, end_adapter, both_adapters):
    '''
    Call cutadapt to remove adapters
    '''

    outfile = os.path.join(get_config_key('general', 'output_directory'), 'trimmed_rna.fastq')
    command = [
        get_config_key('cli-tools', 'trim', 'path_to_trim'),
        file_to_trim
    ]

    if start_adapter is not None:
        command.append('-a')
        command.append(start_adapter)

    if end_adapter is not None:
        command.append('-g')
        command.append(end_adapter)

    if both_adapters is not None:
        command.append('-b')
        command.append(both_adapters)

    command = command + get_config_key('cli-tools', 'trim', 'trim_params')

    print('====> Removing adpters from ends')
    print(command)
    with open(outfile, 'w') as f:
        result = run(command, stdout=f)

    result.check_returncode()

    return result.returncode == 0