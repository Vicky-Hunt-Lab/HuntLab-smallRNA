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