# Copyright 2022 Vicky Hunt Lab Members
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import os.path

from subprocess import run

from .config import do_log, get_config_key

def run_trim(file_to_trim, start_adapter, end_adapter, both_adapters, quiet=0):
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

    if get_config_key('cli-tools', 'trim', 'trim_pass_threads'):
        threads = get_config_key('general', 'threads')

        command.append('-j')
        command.append(str(threads))

    command = command + get_config_key('cli-tools', 'trim', 'trim_params')

    do_log(quiet, '====> Removing adpters from ends')
    do_log(quiet, command)
    with open(outfile, 'w') as f:
        result = run(command, stdout=f)

    result.check_returncode()

    return result.returncode == 0